

rgbImage = imread("special_mat.png");
% Get the dimensions of the image.  numberOfColorBands should be = 3.
[rows, columns, numberOfColorBands] = size(rgbImage);
% Display the original color image.
subplot(2, 1, 1);
imshow(rgbImage);
title('Original Color Image');

hsvImage = rgb2hsv(rgbImage);
hImage = hsvImage(:, :, 1);
sImage = hsvImage(:, :, 2);
vImage = hsvImage(:, :, 3);

mask = vImage == 0; % or whatever.

% Extract the individual red, green, and blue color channels.
redChannel = rgbImage(:, :, 1);
redChannel = regionfill(redChannel, mask);
greenChannel = rgbImage(:, :, 2);
greenChannel = regionfill(greenChannel, mask);
blueChannel = rgbImage(:, :, 3);
blueChannel = regionfill(blueChannel, mask);
% Construct Gaussian Kernel
m=41; 
n=41;
% sigma=15;
% [h1, h2] = meshgrid(-(m-1)/2:(m-1)/2, -(n-1)/2:(n-1)/2);
% h1 = pi*h1./max(h1(:)) - h1./max(h1(:));
% h2 = pi*h2./max(h2(:)) - h2./max(h2(:));
% hg = cos(h1).^2 + cos(h2).^2;
% %hg= exp(-(h1.^2+h2.^2)/(2*sigma^2));            %Gaussian function
% h=hg ./sum(hg(:));
h = ones(m, n)./(m*n);
% % Could be done easier with fspecial though!
% Convolve the three separate color channels.
redBlurred = conv2(redChannel, h, 'same');
greenBlurred = conv2(greenChannel, h, 'same');
blueBlurred = conv2(blueChannel, h, 'same');
% Recombine separate color channels into a single, true color RGB image.
rgbImage2 = cat(3, uint8(redBlurred), uint8(greenBlurred), uint8(blueBlurred));
% Display the blurred color image.
subplot(2, 1, 2);
imshow(rgbImage2);
title('Blurred Color Image');
%imwrite(rgbImage2, "wall_convolved"


