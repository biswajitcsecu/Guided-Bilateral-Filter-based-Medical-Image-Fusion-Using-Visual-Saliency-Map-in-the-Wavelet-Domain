
clc; close all; clear; warning off; format compact

inImg = im2double(rgb2gray(imread('src.jpg')));
inImg = imresize(inImg, [128, 128], 'bilinear');
subplot(1,2,1);
imshow(inImg);title('Original Image');

FFT = fft2(inImg); 
LogAmplitude = log(abs(FFT));
Phase = angle(FFT);
SpectralResidual = LogAmplitude - imfilter(LogAmplitude, fspecial('average', 3), 'replicate'); 
saliencyMap = abs(ifft2(exp(SpectralResidual + 1i*Phase))).^2;

saliencyMap = mat2gray(imfilter(saliencyMap, fspecial('disk', 3)));
subplot(1,2,2);
imshow(saliencyMap, []);title('Saliency Map');
drawnow;



