%Reading our input image
function haze_removal(inputimage)
X=imread(inputimage);
figure,imshow(X);
title('Input Image');
%Our configurable parameters
patch_size=5;
patch_size_1=5;
t0=0.05;
omega=0.95;

%Doing the dark channel extraction
M=min(X,[],3);
rows=size(M,1);
cols=size(M,2);
dark_channel=M;
disp(size(M));
for i=1:1:rows
    for j=1:1:cols
        %place_holder=M(max(1,i-patch_size):min(rows,i+patch_size),max(1,j-patch_size):min(cols,j+patch_size));
        dark_channel(i,j)=min(M(max(1,i-patch_size):min(rows,i+patch_size),max(1,j-patch_size):min(cols,j+patch_size)),[],'all');
    end
end

%Let us do another dark channel extraction here, this is just for trying

M=min(X,[],3);
rows=size(M,1);
cols=size(M,2);
dark_channel_1=M;
disp(size(M));
for i=1:1:rows
    for j=1:1:cols
        %place_holder=M(max(1,i-patch_size):min(rows,i+patch_size),max(1,j-patch_size):min(cols,j+patch_size));
        dark_channel_1(i,j)=min(M(max(1,i-patch_size_1):min(rows,i+patch_size_1),max(1,j-patch_size_1):min(cols,j+patch_size_1)),[],'all');
    end
end

dark_channel_1=imbilatfilt(dark_channel_1,50,40);
%disp(size(dark_channel));

figure,imshow(dark_channel_1);
title('dark channel image');

%Detecting the atmospheric light. Typically the atmospheric light are the
%brightest pixels in the dark channel
disp(find(dark_channel_1==214));

num_brightpixel=uint16(rows*cols/1000);
[check, sortIndex] = sort(dark_channel_1(:), 'descend');  % Sort the values in descending order
maxIndex = sortIndex(1:num_brightpixel);
atmos=[1,1,1];
for i=1:1:num_brightpixel
    temp2=uint16(maxIndex(i)/rows+1);
    temp1=uint16(mod(maxIndex(i),rows));
    %disp(dark_channel_1(temp1,temp2));
    if temp2 == 0 
        temp2=cols;
    end
    if temp1 ==0
        temp1= rows;
    end
    atmos(1)=max(atmos(1),X(temp1,temp2,1));
    atmos(2)=max(atmos(2),X(temp1,temp2,2));
    atmos(3)=max(atmos(3),X(temp1,temp2,3));
    %X(temp1,temp2,:)=[0,0,255];
end
disp(atmos(:));

%Hence we have obtained the atmospheric light
%Now let us do the haze removal
output=X;
transmap=zeros(size(output,1),size(output,2));

%Finding the transmission map
for i=1:1:rows
    for j=1:1:cols
        t1=1-omega*double(dark_channel_1(i,j))/double(atmos(1));
        t2=1-omega*double(dark_channel_1(i,j))/double(atmos(2));
        t3=1-omega*double(dark_channel_1(i,j))/double(atmos(3));
        transmap(i,j)=min([t1 t2 t3]);
    end   
end

%Doing soft matting by billateral filtering
%transmap=imbilatfilt(transmap,50,10);
transmap1 = soft_matting(output, transmap);
%transmap1 = transmap;
%FInding the final output image
for i=1:1:rows
    for j=1:1:cols
        %output(i,j,1)=(double(X(i,j,1))-double(atmos(1)))/double(max(t1,t0))+double(atmos(1));
        output(i,j,1)=(double(X(i,j,1))-double(atmos(1)))/double(max(transmap1(i,j),t0))+double(atmos(1));
        %output(i,j,2)=(double(X(i,j,2))-double(atmos(2)))/double(max(t2,t0))+double(atmos(2));
        output(i,j,2)=(double(X(i,j,2))-double(atmos(2)))/double(max(transmap1(i,j),t0))+double(atmos(2));
        %output(i,j,3)=(double(X(i,j,3))-double(atmos(3)))/double(max(t3,t0))+double(atmos(3));  
        output(i,j,3)=(double(X(i,j,3))-double(atmos(3)))/double(max(transmap1(i,j),t0))+double(atmos(3));  
    end   
end

figure,imshow(output);
title('De-hazed Image');
figure,imshow(transmap);
title('transmap');
end
