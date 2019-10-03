
clc; close all; clear; warning off; format compact

myimage=imread('src.jpg');
image = myimage;
wavename = 'haar';
[cA,cH,cV,cD] = dwt2(im2double(image),wavename);
[cAA,cAH,cAV,cAD] = dwt2(cA,wavename); % Recompute Wavelet of Approximation Coefs.
Level2=[cAA,cAH; cAV,cAD]; %contacinat
figure('Position',[10 10 740 640]);
imshow([Level2,cH; cV,cD],'Colormap',gray); %2 level
drawnow;



