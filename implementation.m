clc;
tic;
clear all;
close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% X-ray image segmentation 
I = imread('Abdul.jpg');
figure
imshow(I)

I = I(:,:,1);

[yl xl] = size(I);
xi = [0, xl];
yi = [yl*0.60, yl*0.60];
[cx, cy, c] = improfile(I,xi, yi); 
flag1 = 1;
for i=1:xl
    if c(i) > 60 && flag1
        x1 = cx(i);
        flag1 = 0;
    end
    if c(i) < 20 && flag1 == 0
        x2 = cx(i);
        break;
    end
end
xi = [0, xl];
yi2 = [yl*0.98, yl*0.98];
[cx, cy, c] = improfile(I,xi, yi2); 

flag1 = 1;
up = 0;
down = 0;
for i=1:yl
    if c(i) < 130
        y2 = cy(i);
        break;
    end
end


I2 = imcrop(I, [x1, yi(1), x2-x1, y2-yi(1)-4]);
figure(2);
imshow(I2)
title('Cropped focused image');
[~, threshold] = edge(I2, 'sobel');
fudgeFactor = .5;
BWs = edge(I2,'sobel', threshold * fudgeFactor);
figure(3)
imshow(BWs)
title('binary gradient mask');
imwrite(BWs, 'exp.jpg');
se90 = strel('line', 2, 90);
se0 = strel('line', 2, 0);
BWsdil = imdilate(BWs, [se90 se0]);
figure(4)
imshow(BWsdil)
title('dilated gradient mask');

sedisk = strel('disk', 4);
bw = imfill(BWsdil, 'holes');
bw2 = ~bwareaopen(~bw, 10);
imshow(bw2)
D = -bwdist(~bw);
imshow(D,[])
Ld = watershed(D);
imshow(label2rgb(Ld))
bw2 = bw;
bw2(Ld == 0) = 0;
imshow(bw2)
mask = imextendedmin(D,2);
imshowpair(bw,mask,'blend')
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = bw;
bw3(Ld2 == 0) = 0;
imshow(bw3)
BWdfill = bw3;
figure(5)
imshow(BWdfill);
title('binary image with filled holes');
BWnobord = imclearborder(BWdfill, 4);
figure(6)
imshow(BWnobord)
title('cleared border image');

radius = 5;
decomposition = 0;
se = strel('disk', radius, decomposition);
BWfinal = imopen(BWnobord, se);
figure(7)
imshow(BWfinal) 
title('final segmented image');

%{
f1=I2;
f=f1(:,:,1);
f=double(adapthisteq(uint8(f)));
thrsld = graythresh(f/255);
I =255*double(imadjust(f/255,[thrsld 1],[0 1]));
figure(8)
imshow(uint8(I));
I = uint8(I); 
%}

%{
figure(9)
imshow(bw)
bw2 = ~bwareaopen(~bw, 10);
imshow(bw2)
D = -bwdist(~bw);
imshow(D,[])
Ld = watershed(D);
imshow(label2rgb(Ld))
bw2 = bw;
bw2(Ld == 0) = 0;
imshow(bw2)
mask = imextendedmin(D,2);
imshowpair(bw,mask,'blend')
D2 = imimposemin(D,mask);
Ld2 = watershed(D2);
bw3 = bw;
bw3(Ld2 == 0) = 0;
figure(10);
imshow(bw3)
title('Final segmented image')
%}


