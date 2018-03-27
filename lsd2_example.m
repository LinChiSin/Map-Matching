% lsd2_example.m
% Test LSD algorithm with MATLAB
%% show the image.
im = imread('F8ÅäÉ«1.jpg');
imshow(im);
%% show the binary image after the process of LSD.
% note: input parameter is the path of image, use '/' as file separator.
figure;
F=lsd2('F8ÅäÉ«1.jpg');
imshow(F);
