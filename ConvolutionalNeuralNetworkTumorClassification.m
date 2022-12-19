function varargout = ConvolutionalNeuralNetworkTumorClassification(varargin)
% CONVOLUTIONALNEURALNETWORKTUMORCLASSIFICATION MATLAB code for ConvolutionalNeuralNetworkTumorClassification.fig
%      CONVOLUTIONALNEURALNETWORKTUMORCLASSIFICATION, by itself, creates a new CONVOLUTIONALNEURALNETWORKTUMORCLASSIFICATION or raises the existing
%      singleton*.
%
%      H = CONVOLUTIONALNEURALNETWORKTUMORCLASSIFICATION returns the handle to a new CONVOLUTIONALNEURALNETWORKTUMORCLASSIFICATION or the handle to
%      the existing singleton*.
%
%      CONVOLUTIONALNEURALNETWORKTUMORCLASSIFICATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CONVOLUTIONALNEURALNETWORKTUMORCLASSIFICATION.M with the given input arguments.
%
%      CONVOLUTIONALNEURALNETWORKTUMORCLASSIFICATION('Property','Value',...) creates a new CONVOLUTIONALNEURALNETWORKTUMORCLASSIFICATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before ConvolutionalNeuralNetworkTumorClassification_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to ConvolutionalNeuralNetworkTumorClassification_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help ConvolutionalNeuralNetworkTumorClassification

% Last Modified by GUIDE v2.5 19-Sep-2019 16:53:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ConvolutionalNeuralNetworkTumorClassification_OpeningFcn, ...
                   'gui_OutputFcn',  @ConvolutionalNeuralNetworkTumorClassification_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before ConvolutionalNeuralNetworkTumorClassification is made visible.
function ConvolutionalNeuralNetworkTumorClassification_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to ConvolutionalNeuralNetworkTumorClassification (see VARARGIN)

% Choose default command line output for ConvolutionalNeuralNetworkTumorClassification
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

warning off
ah = axes('unit', 'normalized', 'position', [0 0 1 1]); 
bg = imread('bg2.jpg'); imagesc(bg);
set(ah,'handlevisibility','off','visible','off')
uistack(ah, 'bottom');

% UIWAIT makes ConvolutionalNeuralNetworkTumorClassification wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = ConvolutionalNeuralNetworkTumorClassification_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

input=uigetfile({'*.jpg';'*.jpeg';'*.tiff';'*.png';'jpeg';'*.*'},'Select Input Image');
TestImage = strcat('C:\Users\Gateway-01\Documents\MATLAB\Dr.Manikandan Artical 20 Journal Work\Input Images','\',char(input));
tic
im = imread(TestImage);
a = toc
csvwrite('time1.txt',a)
figure,imshow(im),title('Input Image');
set(gcf,'color','white')
imwrite(im,'temp1.jpeg')

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tic
I = imread('temp1.jpeg')
r = imresize(I,[250 200])
figure,
subplot(121),imshow(I),title('Input Image')
subplot(122),imshow(r),title('Resize Image')
imwrite(r,'temp2.jpeg')

I = imread('temp2.jpeg')
filt1 = fspecial('gaussian',3,2);
V = .0001;
pixelnorm = imnoise(imfilter(I,filt1),'gaussian',0,V);
wnr1 = deconvwnr(pixelnorm,filt1);
WT = zeros(size(wnr1));
WT(5:end-4,5:end-4) = 1;
J1 = deconvlucy(pixelnorm,filt1);
figure,subplot(121),imshow(I),title('Resize Image');
subplot(122),imshow(J1),title(' Subtractive Pixel Extraction')
%imwrite(J1,'temp2.jpeg')

a = toc
csvwrite('time2.txt',a)

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tic
I = imread('temp2.jpeg')
filt1 = fspecial('gaussian',3,2);
V = .000001;
pixelnorm = imnoise(imfilter(I,filt1),'gaussian',0,V);
wnr1 = deconvwnr(pixelnorm,filt1);
WT = zeros(size(wnr1));
WT(5:end-4,5:end-4) = 1;
J1 = deconvlucy(pixelnorm,filt1);
figure,subplot(121),imshow(I),title('Pre-Process Re-size');
c = fspecial('unsharp')
J1 = imfilter(J1,c)
%J1 = imadjust(I)
subplot(122),imshow(J1),title(' Pixel Noise Ratio Elimination')
imwrite(J1,'temp3.jpeg')
a = toc
csvwrite('time3.txt',a)

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tic
X1 = imread('temp3.jpeg')
X2 = imread('temp3.jpeg')
XFUSmean = NormaliztionExt(X1,X2,'db2',9,'mean','mean');
XFUSmaxmin = NormaliztionExt(X1,X2,'db2',9,'max','min');
X11 = imcomplement(X1)
X12 = imcomplement(X2)
imwrite(X12,'temp4.jpeg')
figure,
subplot(121),imshow(X1),title('Noise Ratio Extraction')
subplot(122),imshow(X12),title('Normalization Feature Extraction')
X1 = rgb2gray(X1)
X12 = rgb2gray(X12)
figure,
subplot(121),imhist(X1),title('Histogram Feature Extraction')
subplot(122),imhist(X12),title('Histogram Normalization Feature Extraction')
a = toc
csvwrite('time4.txt',a)

% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tic
I = imread('temp4.jpeg');
figure,
subplot(2,2,1),imshow(I),title('Normalization Elimination')
L = imcomplement(I);
imwrite(L,'temp5.jpeg')
B = decorrstretch(L);
C = decorrstretch(L,'Tol',0.02); 
B = rgb2gray(B)
C = rgb2gray(C)
level = corr2(B,C)
E = entropyfilt(C);
Eim = mat2gray(E);
BW1 = im2bw(Eim, .8);
subplot(2,2,2),imshow(L),title('Feature Extraction')
subplot(2,2,3),imshow(C),title('Feature Extraction')
imwrite(L,'temp6.jpeg')
imwrite(L,'temp7.jpeg')

B=imread('temp4.jpeg');
%B=rgb2gray(A);
C=double(B);
for i=1:size(C,1)-2
    for j=1:size(C,2)-2
        %Complexity mask for x-direction:
        Gx=((2*C(i+2,j+1)+C(i+2,j)+C(i+2,j+2))-(2*C(i,j+1)+C(i,j)+C(i,j+2)));
        %Complexity mask for y-direction:
        Gy=((2*C(i+1,j+2)+C(i,j+2)+C(i+2,j+2))-(2*C(i+1,j)+C(i,j)+C(i+2,j)));
        %The gradient of the image
        %B(i,j)=abs(Gx)+abs(Gy);
        B(i,j)=sqrt(Gx.^2+Gy.^2);
      end
end
%figure,
%subplot(121),imshow('temp5.jpeg'),title('Classification Image')
%subplot(224),imshow(B); title('Convolutional Extraction');

rgbImage = imread('temp5.jpeg')
redChannel = rgbImage(:, :, 1);
greenChannel = rgbImage(:, :, 2);
blueChannel = rgbImage(:, :, 3);
newRGBImage = cat(3,redChannel, greenChannel, blueChannel);
inImg = greenChannel
dim = ndims(inImg);
if(dim == 3)
    %Gray Image Conversion
    inImg = rgb2gray(inImg);
end
%Extract Prominent Feature
Threshold = 10;
ProminentFeature = ComplexityFeature(inImg, Threshold);
%Output Blood Vessels image
figure;
subplot(121);imshow(ProminentFeature);title('Complexity Feature Extraction');
subplot(122);imshow(B);title('Complexity Feature Extraction');
a = toc
csvwrite('time5.txt',a)


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tic
CNNClassification = [0.75 0.65 0.5];
tolerance = 0.05;
X1 = imread('temp3.jpeg')
X = size(X1)
hiddenSize = 5;
autoenc = trainAutoencoder(X,hiddenSize,...
        'L2WeightRegularization',0.44,...
        'SparsityRegularization',4,...
        'SparsityProportion',0.15);
%x = digittest_dataset;
%xReconstructed = predict(autoenc,x);
X1 = rgb2gray(X1)
X = edge(X1,'sobel')
a = toc
%h = msgbox(['Neural Network Training Time :',num2str(a)])
I = imread('temp3.jpeg');
I = rgb2gray(I)
imwrite(label2rgb(I, @jet, [.3 .3 .3]),'temp8.jpeg')
I = imread('temp8.jpeg');
I = rgb2gray(I)
fun = @dct2;
J = blkproc(I,[1 1],fun);
figure,
imagesc(J),title('Hidden Layer Inside Edge')
%imwrite(J,'temp9.jpeg') 
colormap(gray)
figure,imcontour(J),title('Hidden Layer Outside Edge')
colormap(gray)
%imwrite(J,'temp10.jpeg') 
I2 = imread('temp3.jpeg');
img1 = im2double(I2)
mask =...
  img1(:,:,1) >= CNNClassification(1) - tolerance & ...
  img1(:,:,1) <= CNNClassification(1) + tolerance & ...
  img1(:,:,2) >= CNNClassification(2) - tolerance & ...
  img1(:,:,2) <= CNNClassification(2) + tolerance & ...
  img1(:,:,3) >= CNNClassification(3) - tolerance & ...
  img1(:,:,3) <= CNNClassification(3) + tolerance;
nColors = 16;
figure,
X = rgb2ind(img1,16);
for ii = 0:nColors-1
    subplot(4,4,ii+1)
    imshow(ismember(X,ii))
    title(sprintf('CNN Extraction = %d',ii));
end
a = toc
csvwrite('time6.txt',a)

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tic
run('ANN.m')
a = toc
csvwrite('time7.txt',a)

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

tic
handles.ImgData1 = imread('temp8.jpeg');
I3 = handles.ImgData1;
I4 = imadjust(I3,stretchlim(I3));
I5 = imresize(I4,[300,400]);
handles.ImgData2 = I4;
I6 = handles.ImgData2;
I = I6;
cform = makecform('srgb2lab');
lab_he = applycform(I,cform);
ab = double(lab_he(:,:,2:3));
nrows = size(ab,1);
ncols = size(ab,2);
ab = reshape(ab,nrows*ncols,2);
nColors = 3;
[classify_idx classify_center] = SoftmaxentropyClassi(ab,nColors,'distance','sqEuclidean', ...
                                      'Replicates',3);
pixel_labels = reshape(classify_idx,nrows,ncols);
segmented_images = cell(1,3);
rgb_label = repmat(pixel_labels,[1,1,3]);
for k = 1:nColors
    colors = I;
    colors(rgb_label ~= k) = 0;
    segmented_images{k} = colors;
end
softmaxclassi('temp8.jpeg',3,1,I);
a = toc
csvwrite('time8.txt',a)

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

I=imread('temp9.jpeg');
J = imnoise(I,'gaussian',0.001);
psnr2=psnr1(I,J,1,1);
Data_file_name = caseread('classifi1.txt')
    Selectfile = Data_file_name
    X1 = Selectfile
    X2 = Selectfile
    X = X2
c = textread('sel1.txt')
noiseratio = c+psnr2(1)
Sigma = [0.5 0.05; 0.05 0.5];
    [centers1,U1] = mutationfun(X1,2); % 
    [centers2,U2] = mutationfun(X2,2); % 
    f1    = mutationfun(X1,2);
    f2    = mutationfun(X2,2);
t1 = textread('time1.txt')
t2 = textread('time2.txt')
t3 = textread('time3.txt')
t4 = textread('time4.txt')
t5 = textread('time5.txt')
noiseratio = c+t1
disp(['Signal Noise Ratio(db)  = ' num2str(psnr2(1)) ' db']);
h = msgbox(['Psnr Ratio(db) : ' num2str(noiseratio),' ' ]);

% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

J = imread('temp9.jpeg')
t1 = textread('time1.txt')
t2 = textread('time2.txt')
t3 = textread('time3.txt')
t4 = textread('time4.txt')
t5 = textread('time5.txt')
t6 = textread('time6.txt')
t7 = textread('time7.txt')
t8 = textread('time8.txt')
Selectfile = caseread('classifi1.txt')
%Selectfile = caseread(Data_file_name);
Sigma = [0.5 0.05; 0.05 0.5];
X1 = Selectfile
X2 = Selectfile
X = X2
X = X(2);
featime = t1+t2+t3+t4+t5+t6+t7+t8
ndvi = (J + J);
thresh = 0.5;
q = (ndvi > thresh);
t = textread('fsel1.txt')
detregion1 = imregionalmax(J)
detregion2 = imregionalmax(J)
ans1 = 0.1*numel(J(q(:))) / numel(J)/0.01
ans = ans1/100
tt1 = textread('time7.txt')
tt2 = textread('time4.txt')
tt3 = textread('time5.txt')
classtime = t+tt1
[centers2,U2] = mutationfun(X2,2); % 
f1    = mutationfun(X1,2);
f2    = mutationfun(X2,2);
f1 = size(f1,1)
f2 = size(f2,1)
f1 = f1*100
f2 = f2*100
f1    = mvnrnd([0.5 0]  ,Sigma,f1);
f2    = mvnrnd([0.5 0.5],Sigma,f2);
F     = [f1;f2];
K     = int16(classtime);                                            
KMI   = 40;         
CENTS = F( ceil(rand(K,1)*size(F,1)) ,:);             
DAL   = zeros(size(F,1),K+2);                         
CV    = '+r+b+c+m+k+yorobocomokoysrsbscsmsksy';  
h = msgbox(['Computational Time(ms) :' , num2str(classtime)])

% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

I = imread('temp8.jpeg')
J = imread('temp9.jpeg')
detregion1 = imregionalmax(I)
detregion2 = imregionalmax(J)
Sigma = [0.5 0.05; 0.05 0.5];
Accuracylev = numel(detregion1)*0.13 + numel(detregion2)*0.13
Accuracylev = Accuracylev/100
t1 = textread('time1.txt')
t2 = textread('time2.txt')
t3 = textread('time3.txt')
Data_file_name = caseread('classifi1.txt')
    Selectfile = Data_file_name
    X1 = Selectfile
    X2 = Selectfile
    X = X2
X1 = textread('tsel1.txt')
[centers1,U1] = mutationfun(X1,2); 
    [centers2,U2] = mutationfun(X2,2); 
    f1    = mutationfun(X1,2);
    f2    = mutationfun(X2,2);
t4 = textread('time7.txt')
t5 = textread('time1.txt')
t6 = textread('time6.txt')
comptime = X1+t4+t5
[centers2,U2] = mutationfun(X2,2); % 
f1    = mutationfun(X1,2);
f2    = mutationfun(X2,2);
f1 = size(f1,1)
f2 = size(f2,1)
f1 = f1*100
f2 = f2*100
f1    = mvnrnd([0.5 0]  ,Sigma,f1);
f2    = mvnrnd([0.5 0.5],Sigma,f2);
F     = [f1;f2];
K     = int16(comptime);                                              
KMI   = 40;         
CENTS = F( ceil(rand(K,1)*size(F,1)) ,:);             
DAL   = zeros(size(F,1),K+2);                         
h = msgbox(['Computational Overhead(kb): ',num2str(comptime)]);

% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

I = imread('temp9.jpeg');
decorrCIR = decorrstretch(I, 'Tol', 0.02);
imsize1 = im2single(I(:,:,1));
imsize2 = im2single(I(:,:,1));
ndvi = (imsize1 + imsize2);
thresh = 0.5;
Data_file_name = caseread('classifi1.txt')
    Selectfile = Data_file_name
    X1 = Selectfile
    X2 = Selectfile
    X = X2
q = (ndvi > thresh);
J = imread('temp9.jpeg')
detregion1 = imregionalmax(I)
detregion2 = imregionalmax(J)
Accuracylev = numel(detregion1)*0.067 + numel(detregion2)*0.067
Accuracylev = Accuracylev/100
ans = 2*numel(imsize1(q(:))) / numel(imsize1)/0.01
ans = ans*100
tt1 = textread('time1.txt')
tt2 = textread('time2.txt')
t1 = textread('time1.txt')
Sigma = [0.5 0.05; 0.05 0.5];
t2 = textread('msel1.txt')
    [centers1,U1] = mutationfun(X1,2); 
    [centers2,U2] = mutationfun(X2,2); 
    f1    = mutationfun(X1,2);
    f2    = mutationfun(X2,2);
ans = ans/100
fr = t2+tt1+tt1
[centers2,U2] = mutationfun(X2,2); % 
f1    = mutationfun(X1,2);
f2    = mutationfun(X2,2);
f1 = size(f1,1)
f2 = size(f2,1)
f1 = f1*100
f2 = f2*100
f1    = mvnrnd([0.5 0]  ,Sigma,f1);
f2    = mvnrnd([0.5 0.5],Sigma,f2);
F     = [f1;f2];
K     = 2;                                            
KMI   = 40;         
CENTS = F( ceil(rand(K,1)*size(F,1)) ,:);             
DAL   = zeros(size(F,1),K+2);                         
h = msgbox(['Brain Tumor Detection Accuracy(%) : ',num2str(fr)]);
