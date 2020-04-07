function [ matted_transmap ] = soft_matting( inputimage, transmap )
%bar = waitbar(0,'Matting ..');
rows = size(inputimage,1);
cols = size(inputimage,2);
img_size = rows*cols;
%imgcomp = size(inputimage,3); 
lambda = 1e-4;
epsilon = 1e-8;
U3 = eye(3);
transmap_vector = reshape(transpose(transmap), img_size, 1);
L = sparse(img_size, img_size);
L_row = zeros((img_size + cols)*9, 1); %each element of window stored separately, summed at the end by sparsing!!
L_col = L_row;
L_values = L_col; % function value
a = 1; %index value
for i = 2:1:rows-1
 for j = 2:1:cols-1
  windowimage=inputimage(i-1:i+1, j-1:j+1 , :);
  windowimage_vector = reshape(windowimage, 9, 3); % size = 9 x 3
  %windowsize = size(windowimage, 1)*size(windowimage, 2);
  %windowsize = size(window_index_vector, 1);
  %disp(windowsize);
  window_mean = mean(windowimage_vector); % 3 X 1
  %disp(size(window_mean));
  %mean_windowimage_vector = double(windowimage_vector) - window_mean;
  window_covariance = cov(double(windowimage_vector), 1);
  inverse = window_covariance + epsilon.*U3./9; % size 3 x 3
  %disp(size(inverse));
  for k = i-1:1:i+1
      for p = j-1:1:j+1  %mn iterations for each Ii image pixel
       L_row_value = (k-1)*cols + p; % image index from 2 == window index 1
       for m = i-1:1:i+1
          for n = j-1:1:j+1 %mn iterations for each Ij image pixel
            L_col_value = (m-1)*cols + n; % same as above
            windowimage_vector = reshape(inputimage(k,p,:), 3, 1); % size = 3 x 1
            %disp(size(windowimage_vector));
            mean_windowimage_vector_pixel_i = double(windowimage_vector) - transpose(window_mean); % Same notation as in the Paper
            %disp(size(mean_windowimage_vector_pixel_i));
            windowimage_vector = reshape(inputimage(m,n,:), 3, 1); % size = 3 x 1
            mean_windowimage_vector_pixel_j = double(windowimage_vector) - transpose(window_mean);
            L_row(a) = L_row_value;
            L_col(a) = L_col_value;
            L_values(a) = (eq(L_row_value, L_col_value) - (1 + (transpose(mean_windowimage_vector_pixel_i) / inverse )* mean_windowimage_vector_pixel_j)./9 );
            a = a + 1; %matlab performs division efficiently than inverse!!
          end
       end
      end
  end
 end
%disp(imageRow); 
%waitbar(imageRow/imageSize(1), bar);
end
L = sparse(L_row, L_col, L_values, img_size, img_size); %Summation of same indices here directly performed here!
%(L+xU)t = xt~ => (L+xU)t(tT) = (L+xU) = xt~tT => tT = (L+xU)/xt
matted_transmap_transpose = (L + lambda*speye(img_size))\(lambda*transmap_vector); %Matrix left division!
matted_transmap = transpose(reshape(matted_transmap_transpose, cols, rows)); %tramspose, so reshape it to its transpose!!
end
