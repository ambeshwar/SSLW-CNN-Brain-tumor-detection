I = imread('temp1.jpeg')
v = fspecial('unsharp')
v = imadjust(stretchlim(v))
d = imfilter(I,v)
%d1 = imresize(d,[360 640])
rs = imadjust(I,stretchlim(I));
figure
subplot(121),imshow(I),title('Quality Improve Image');
subplot(122),imshow(rs),title('Quality Improve Image');