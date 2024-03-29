clear all
close all

num_lin = 500;
num_ang = 500;
margin = 1;
mirror_BRDF = @(angle_diff, margin) (1*(abs(angle_diff) <= margin/2));
sigma = 1;
blurred_mirror_BRDF = @(angle_diff, sigma) (normpdf(-angle_diff, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 15, sigma)/normpdf(0, 0, sigma) + normpdf(-angle_diff, -35, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 35, sigma)/normpdf(0, 0, sigma));
blurred_mirror_BRDFInterf = @(angle_diff, sigma) (xcorr(normpdf(-angle_diff, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 15, sigma)/normpdf(0, 0, sigma) + normpdf(-angle_diff, -35, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 35, sigma)/normpdf(0, 0, sigma)));


obs_pos = [0, -5];
obs_size = 1;
obs_interval = [-45, 45];
obs = build_obs(obs_pos, obs_size, obs_interval);

source_pos = [0, -5];
source_support_width = 10;
source_support_size = 2;
[source, interf_test] = build_source(source_pos, source_support_width, source_support_size, obs_size, false);
[empty_source, ~] = build_source(source_pos, source_support_width, source_support_size, obs_size, true);

%% Rendering calls
% Reference rendering
% tic
[x0, h0, g0] = FullRendering(obs, source, blurred_mirror_BRDF, num_lin, num_ang, margin);
[x1, h1, g1] = FullRendering(obs, source, blurred_mirror_BRDF, num_lin, num_ang, sigma);
% toc
% g0 = g0';
% % Single convolution approximation
% tic
g2 = FastRendering(obs, reshape(source(:,3,:), source_support_size, obs_size), h0, num_ang, empty_source);
% toc
% g1 = g1';
% % Rendering operator
% tic
% g2 = FastRendering(obs, reshape(source(:,3,:), source_support_size, obs_size), h0, num_ang, empty_source);
% toc
% % Interferometric rendering from interferometric quantities
% brdf_interf = xcorr(h0);
% [temp1, indexes] = PrepFastRenderingInterf(obs, num_ang, empty_source, source_support_size, num_lin);
% tic
% g3 = FastRenderingInterf(interf_test, brdf_interf, num_ang, obs_size*obs_size, num_lin, indexes, temp1); 
% toc
% % Rendering using matrix product
% g4 = zeros(num_lin, obs_size);
% tic
% mat = CreateRenderingMatrixFromBRDF(obs, source, blurred_mirror_BRDF, num_lin, sigma);
% mat = mat';
% for i=1:obs_size    
%     g4(:, i) = mat*source(:, 3, i);
% end
% toc

% figure;
% subplot(1, 3, 1)
% stem(source(:, 3, 1))
% title("Light source");
% ylabel("Intensity");
% xlabel("m");
% subplot(1, 3, 2)
% stem(linspace(-180, 180, num_ang), h0)
% title("BRDF");
% ylabel("\rho");
% xlabel("Angle difference");
% subplot(1, 3, 3)
% stem(x0, g0)
% title("Rendered signal");
% ylabel("Reflected intensity");
% xlabel("Wall points (x_w)");

figure;
subplot(1, 2, 1)
stem(x1, g1)
title("Rendered signal");
ylabel("Reflected intensity");
xlabel("Wall points (x_w)");
subplot(1, 2, 2)
stem(x1, g2)
title("Rendered signal (Spherical conv.)");
ylabel("Reflected intensity");
xlabel("Wall points (x_w)");
% %% Plots
% figure;
% subplot(1, 2, 1)
% imagesc(g0);
% subplot(1, 2, 2)
% imagesc(g1);
% title("Convolution approx")
% 
% figure;
% subplot(1, 2, 1)
% imagesc(g0);
% subplot(1, 2, 2)
% imagesc(g2);
% title("Full convolution")
% 
% figure;
% subplot(1, 2, 1)
% imagesc(xcorr(g0));
% subplot(1, 2, 2)
% imagesc(g3);
% title("Interf from interf")
% 
% figure;
% subplot(1, 2, 1)
% imagesc(g0);
% subplot(1, 2, 2)
% imagesc(g4);
% title("Matrix product")
% 
% figure;
% subplot(1, 2, 1)
% imagesc(test);
% subplot(1, 2, 2)
% imagesc(g5);
% title("Matrix interf")

%% Functions

function [source, interf_test] = build_source(source_pos, source_support_width, source_support_size, obs_size, isempty)
    source = zeros(source_support_size, 3, obs_size);
    interf_test_prep = zeros(source_support_size, obs_size);
    if(source_support_size == 1)
        x_axis = zeros(1, 1);
    else
        x_axis = linspace(-source_support_width/2, source_support_width/2, source_support_size);
    end
    for i=1:obs_size
        for j=1:source_support_size
            if isempty
                val = 1;
            else
                val = rand();
            end            
            source(j, :, i) = [source_pos(1)+x_axis(j), source_pos(2), val];
            interf_test_prep(j, i) = val;
        end
        if ~isempty
            source(:, 3, i) = 0.9*source(:, 3, i);
            ind = randi(source_support_size);
            source(ind, 3, i) = 1;
        end
    end
    interf_test = xcorr(interf_test_prep);
end

function obs = build_obs(obs_pos, obs_size, obs_interval)
    obs = [obs_pos(1) obs_pos(2) obs_interval(1) obs_interval(2) obs_size];    
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end
