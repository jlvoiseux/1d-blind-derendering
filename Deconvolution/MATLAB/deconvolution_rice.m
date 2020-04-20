clear all
close all

f = double(imread('cameraman.tif'));
%f = rgb2gray(RGB);
%[g_cropped, g, h_toep] = convolve_prod(f, 2.5, 7);
%f_est = tik_deconv(y, h_toep, 0.0001, size(f, 1), size(f, 2));
[g_cropped, g, h_size] = convolve_standard(f, 2.5, 7);
[f_est, h_est] = blind_dconv(g, 100, h_size(1), 0.001, g_cropped, size(g_cropped, 1));
subplot(2,2,1);
imagesc(f)
subplot(2,2,2);
imagesc(g_cropped)
subplot(2,2,3);
imagesc(f_est)
subplot(2,2,4);
imagesc(h_est)
colormap('gray');

function [blurred_cropped, blurred, psf_size] = convolve_standard(img, sigma, psf_size)
    N = psf_size;
    [x, y] = meshgrid(round(-N/2):round(N/2), round(-N/2):round(N/2));
    psf = exp(-x.^2/(2*sigma^2)-y.^2/(2*sigma^2));
    psf = psf./sum(psf(:));
    imagesc(psf)
    blurred = conv2(img, psf);
    blurred_cropped = conv2(img, psf, 'same');
    psf_size = size(psf);
end

function[deconv,psf]=blind_dconv(z,iterations,hesize,alpha,fe,k)
    if(mod(hesize,2)~=1)
      error('Size of psf must be odd');
    end

    YSIZE=size(z,1);
    XSIZE=size(z,2);
    SIZEMAX=max([XSIZE YSIZE]);

    %z is the linearly degraded image, fe is the image estimate, and k is the
    %finite support. he is the PSF estimate.

    %This is the PSF ESTIMATE (we estimate it as an impulse)
    he=zeros(1,hesize);
    he(ceil(hesize/2))=1;
    he=he'*he;

    G=fft2(z,(SIZEMAX),(SIZEMAX)); %z is the degraded image

    n=0;

    while(n<iterations)
      Fe=fft2(fe,(SIZEMAX),(SIZEMAX));Hoe=fft2(he,(SIZEMAX),(SIZEMAX));

      %impose fourier constraints
      Hne=(G.*conj(Fe))./(abs(Fe).^2+(alpha./(abs(Hoe).^2)));
      hne=ifft2(Hne);

      %impose blur constraints
      hne=abs(hne(1:hesize,1:hesize));
      Hoe=fft2(hne,(SIZEMAX),(SIZEMAX));

      %Impose Fourier constraints
      Fne=(G.*conj(Hoe))./((abs(Hoe).^2+(alpha./abs(Fe).^2)));
      fne=ifft2(Fne);

      %impose image constraints
      fne=abs(fne(1:YSIZE,1:XSIZE));
      fne(~k)=0;
      fe=fne;he=hne;
      n=n+1;
    end

    deconv=fne;
    psf=hne;
end