clear all
close all

f = double(imread('rendered.png'));
f1 = f(:, :, 1);
f2 = f(:, :, 2);
f3 = f(:, :, 3);

%% Standard 
%kernel
imgx = size(f1, 1);
imgy = size(f1, 2);
m=10; 
n=10;
[h1, h2] = meshgrid(-(m-1)/2:(m-1)/2, -(n-1)/2:(n-1)/2);
h1 = pi*h1./max(h1(:)) - h1./max(h1(:));
h2 = pi*h2./max(h2(:)) - h2./max(h2(:));
hg = cos(h1).^2 + cos(h2).^2;
%psf=hg ./sum(hg(:));
psf = ones(m, n)./(m*n);
h = convmtx2(psf, imgx, imgy);
disp("toep")

% deconv
f_est1 = tik_deconv(padarray(f1,[9 9],'replicate','post'), h, 1, size(f1, 1), size(f1, 2));
disp("1")
f_est2 = tik_deconv(padarray(f2,[9 9],'replicate','post'), h, 1, size(f2, 1), size(f2, 2));
disp("2")
f_est3 = tik_deconv(padarray(f3,[9 9],'replicate','post'), h, 1, size(f3, 1), size(f3, 2));
disp("3")
final = cat(3, uint8(f_est1), uint8(f_est2), uint8(f_est3));
imagesc(final);

%% Blind
% [f_est1b, h_est1b] = blind_deconv(f1, f1, 100, 100, size(f1, 1), size(f1, 2), 10, 100);
% [f_est2b, h_est2b] = blind_deconv(f2, f2, 100, 100, size(f2, 1), size(f2, 2), 10, 100);
% [f_est3b, h_est3b] = blind_deconv(f3, f3, 100, 100, size(f3, 1), size(f3, 2), 10, 100);
% finalb = cat(3, uint8(f_est1b), uint8(f_est2b), uint8(f_est3b));
% imagesc(finalb);


%% Blind Deconvolution
%[g_cropped, g, h, h_size] = convolve_standard(f, 2.5, 7);
%[f_est, h_est] = blind_deconv(g, g, h_size(1), h_size(2), size(f, 1), size(f, 2), 10, 100);
%colormap('gray');

%% Standard Deconvolution
% [g_cropped, g, h_toep] = convolve_prod(f, 2.5, 7);
% f_est = tik_deconv(g, h_toep, 1, size(f, 1), size(f, 2));
% colormap('gray');


%% Blind Deconvolution
function [blurred_cropped, blurred, psf, psf_size] = convolve_standard(img, sigma, psf_size)
    N = psf_size;
    [x, y] = meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
    psf = exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
    psf = psf./sum(psf(:));
    imagesc(psf);
    blurred = conv2(img, psf);
    %[psf_h psf_w] = size(psf);
    %[img_h img_w] = size(img);
    %blurred = ifft2(fft2(img, img_h+psf_h-1, img_w+psf_w-1).*fft2(psf, img_h+psf_h-1, img_w+psf_w-1));
    blurred_cropped = conv2(img, psf, 'same');
    psf_size = size(psf);
end

function out = apply_fourier_constraints(A, B, C, alpha)
    out = (A.*conj(B))./(abs(B).^2+(alpha./(abs(C).^2)));
end

function [f_est_final, h_est_final] = blind_deconv(blurred_img, guess_img, psf_h, psf_w, img_h, img_w, alpha, count)
    G = fft2(blurred_img);
    Gh = size(G, 1);
    Gw = size(G, 2);
    f_est = guess_img;
    h_est = zeros(psf_h, psf_w);
    h_est(ceil(size(h_est, 1)/2), ceil(size(h_est, 2)/2)) = 1;
    H_est = fft2(h_est, Gh, Gw);
    for i=1:count
        F_est = fft2(f_est, Gh, Gw);
        H_tilde = apply_fourier_constraints(G, F_est, H_est, alpha);
        h_tilde = ifft2(H_tilde);
        h_est = abs(h_tilde(1:psf_h, 1:psf_w));
        H_est = fft2(h_est, Gh, Gw);
        F_tilde = apply_fourier_constraints(G, H_est, F_est, alpha);
        f_tilde = ifft2(F_tilde);
        f_est = abs(f_tilde(1:img_h, 1:img_w));
    end
    f_est_final = f_est;
    h_est_final = h_est;
end

%% Standard Deconvolution
function [blurred_cropped, blurred, h] = convolve_prod(img, sigma, psf_size)
    N = psf_size;
    [x, y] = meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
    psf = exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
    psf = psf./sum(psf(:));   
    imagesc(psf)
    [blurred, h] = convToeplitz(img,psf);
    blurred_cropped = blurred(floor(psf_size/2)+1:end-ceil(psf_size/2)-1, ceil(psf_size/2)+1:end-floor(psf_size/2)-1);
end

function [res, h] = convToeplitz(img, psf)
    imgx = size(img, 1);
    imgy = size(img, 2);
    h = convmtx2(psf, imgx, imgy); % Image processing toolbox
    res = reshape(h*img(:), size(psf)+[imgx imgy]-1);    
end

function estimate = tik_deconv(blurred_img, h, alpha, orig_y, orig_x)
    temp1 = h'*h + alpha*speye(size(h, 2));
    temp2 = h'*blurred_img(:);
    estimate = temp1\temp2;
    estimate = reshape(estimate, [orig_y orig_x]);
end