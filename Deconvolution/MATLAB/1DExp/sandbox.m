clear all
close all

x = linspace(-1, 1, 1000);
mu = 0;
sigma = 0.01;
h_prep1 = normpdf(x, mu, sigma)/normpdf(0, mu, sigma);
delta_num = 5;
h_prep2 = zeros(1, 1000);
h_prep2(1+length(h_prep2)/delta_num/2:length(h_prep2)/delta_num:end) = 1;
h = conv(h_prep1, h_prep2, "same");


delta_num = 3;
f = zeros(1, 1000);
f(1+length(f)/delta_num/2:length(f)/delta_num:end) = 1;
subplot(1, 2, 1)
stem(f);
title('Original signal')
subplot(1, 2, 2);
stem(h);
title('Impulse response')

figure;
g = conv(h, f);
g = g/sum(g);
stem(g);
title('Convolved signal')

figure;
[f_est, h_est] = LucyRichardsonBlind(g, length(f), length(h), 100, 10);
subplot(1, 2, 1)
stem(f_est);
subplot(1, 2, 2);
stem(h_est);