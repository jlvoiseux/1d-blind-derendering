s = rand(1, 100);
h1 = rand(1, 5);
h2 = rand(1, 5);
u1 = conv(s, h1);
u2 = conv(s, h2);
S = convmtx(xcorr(s), 9);
test1 = xcorr(u1, u2);
test2 = S'*(xcorr(h1,h2) + xcorr(h2,h1))';

figure;
subplot(1, 2, 1);
stem(test1./max(test1));
subplot(1, 2, 2);
stem(test2./max(test2));