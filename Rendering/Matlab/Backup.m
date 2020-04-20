tic
r1 = 2.*rand(3, 2, 2, samples, h, w);
disp("r1 " + toc);
tic
r2 = 2.*rand(3, 2, 2, samples, h, w);
disp("r2 " + toc);
tic
dx = ((r1<1).*(sqrt(r1)-1))+((1-(r1<1)).*(1-sqrt(2-r1)));
disp("dx " + toc);
tic
dy = ((r2<1).*(sqrt(r2)-1))+((1-(r2<1)).*(1-sqrt(2-r2)));
disp("dy " + toc);
tic
clear r1;
clear r2;
sxPrep = [0 1 ; 0 1 ; 0 1];
syPrep = [0 0 ; 0 0 ; 0 0];
sxPrep(:, :, 2) = [0 1 ; 0 1 ; 0 1];
syPrep(:, :, 2) = [1 1 ; 1 1 ; 1 1];
sx = repmat(sxPrep, 1, 1, 1, samples, h, w);
disp("sx " + toc);
tic
sy = repmat(syPrep, 1, 1, 1, samples, h, w);
disp("sy " + toc);
tic
x = permute(repmat(repmat(0:w-1, h, 1)', 1, 1, samples, 2, 2, 3), [6, 5, 4, 3, 2, 1]);
disp("x " + toc);
tic
y = permute(repmat(repmat((0:h-1)', 1, w)', 1, 1, samples, 2, 2, 3), [6, 5, 4, 3, 2, 1]);
disp("y " + toc);
tic
camX = repmat(camXconst', 1, 2, 2, samples, h, w);
disp("camX " + toc);
tic
camY = repmat(camYconst', 1, 2, 2, samples, h, w);
disp("camY " + toc);
tic
camD = repmat(camera(2,:)', 1, 2, 2, samples, h, w);
disp("camD " + toc);
d = camX.*(((sx+0.5+dx)./2 + x)./w - 0.5) + camY.*(((sy+0.5+dy)./2 + y)./h - 0.5) + camD; 
disp("d " + toc);