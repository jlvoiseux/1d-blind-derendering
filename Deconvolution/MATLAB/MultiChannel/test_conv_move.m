clear all
close all

x1 = rand(1, 25);
move_vec = [0 1];
test1 = conv(x1, move_vec);
figure;
subplot(1, 2, 1)
stem(x1);
subplot(1, 2, 2)
stem(test1(1:end-1));
