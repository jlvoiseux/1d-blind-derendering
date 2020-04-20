clear all
close all

wiener = true;
lucy = false;
tikh = false;
blind = false;

f = double(imread('todeconvolvesmall.png'));
f1 = f(:, :, 1);
f2 = f(:, :, 2);
f3 = f(:, :, 3);

m=40; 
n=40;
psf = ones(m, n)./(m*n);
% deconv

if lucy == true
    iter = 10;
    f_est1 = deconvlucy(f1, psf);
    f_est2 = deconvlucy(f2, psf);
    f_est3 = deconvlucy(f3, psf);   
elseif wiener == true % nsr
    nsr = 0.2;
    f_est1 = deconvwnr(f1, psf, nsr);
    f_est2 = deconvwnr(f2, psf, nsr);
    f_est3 = deconvwnr(f3, psf, nsr); 
    f2 = f2/(1+nsr);
elseif tikh == true % np
    nsr = 0.2;
    f_est1 = deconvreg(f1, psf, nsr);
    f_est2 = deconvreg(f2, psf, nsr);
    f_est3 = deconvreg(f3, psf, nsr); 
elseif blind == true
    iter = 10;
    f_est1 = deconvblind(f1, psf, iter);
    f_est2 = deconvblind(f2, psf, iter);
    f_est3 = deconvblind(f3, psf, iter); 
end
initial = cat(3, uint8(f1), uint8(f2), uint8(f3));
final = cat(3, uint8(f_est1), uint8(f2), uint8(f_est3));
imagesc(final);
imwrite(final, "blind.png");