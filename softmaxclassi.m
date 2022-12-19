function softmaxclassi(target,n,threshtype,im)
I = imread('temp3.jpeg'); % image
Level = I>210; 
Level_filled = imfill(Level,'holes'); 
DetectionM = and(Level,Level_filled); 
Detection = uint8(DetectionM).*I; 
Detection=double(Detection);
figure,
ax(2)=subplot(222);
imshow(Detection);
title('Normalization Extraction');
%imwrite(Detection,'temp10.jpeg')
d = imread('temp3.jpeg')
l = imresize(Detection,[250 250])
c = size(l)
subplot(221)
imshow(d),title('Normalization Extraction')
subplot(222)
imshow('temp8.jpeg'),title('Complexity Feature Extraction')
subplot(223)
imshow(l),title('Brain Tumor Detection Area')
csvwrite('size2.txt',c)
k = csvread('size2.txt')

