clear all
close all

x1 = rand(1, 100);
x2 = rand(1, 100);

test1 = xcorr(conv(x1, x2), 'normalized');

test2 = conv(xcorr(x1), xcorr(x2));
test2 = test2./max(test2);

figure;
subplot(1, 2, 1);
stem(test1);
subplot(1, 2, 2);
stem(test2);

err = immse(test1, test2);