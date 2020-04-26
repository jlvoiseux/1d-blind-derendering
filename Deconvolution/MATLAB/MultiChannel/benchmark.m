clear all
close all

num_lin = 500;
num_ang = 500;
margin = 5;
mirror_BRDF = @(angle_diff, margin) (1*(abs(angle_diff) <= margin/2));
sigma = 1;
blurred_mirror_BRDF = @(angle_diff, sigma) (normpdf(-angle_diff, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 15, sigma)/normpdf(0, 0, sigma) + normpdf(-angle_diff, -35, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 35, sigma)/normpdf(0, 0, sigma));

obs_pos = [0, -5];
obs_size = 4;
obs_interval = [-45, 45];
obs = build_obs(obs_pos, obs_size, obs_interval);

source_pos = [0, -50];
source_support_width = 10;
source_support_size = 2;
[source, interf_test] = build_source(source_pos, source_support_width, source_support_size, obs_size, false);
[empty_source, ~] = build_source(source_pos, source_support_width, source_support_size, obs_size, true);

%% Rendering calls
% Reference rendering
tic
[x0, h0, g0] = FullRendering(obs, source, blurred_mirror_BRDF, num_lin, num_ang, sigma);
toc
g0 = g0';
% Single convolution approximation
tic
g1 = RenderingConv(obs, source, blurred_mirror_BRDF, sigma, num_lin, num_ang, mirror_BRDF, margin);
toc
g1 = g1';
% Rendering operator
tic
g2 = FastRendering(obs, reshape(source(:,3,:), source_support_size, obs_size), h0, num_ang, empty_source);
toc
% Interferometric rendering from interferometric quantities
tic
brdf_interf = xcorr(h0);
g3 = FastRenderingInterf(obs, interf_test, brdf_interf, num_ang, empty_source);
toc

%% Plots
figure;
subplot(1, 2, 1)
imagesc(g0);
subplot(1, 2, 2)
imagesc(g1);
title("Convolution approx")

figure;
subplot(1, 2, 1)
imagesc(g0);
subplot(1, 2, 2)
imagesc(g2);
title("Full convolution")

figure;
subplot(1, 2, 1)
imagesc(xcorr(g0));
subplot(1, 2, 2)
imagesc(g3);
title("Interf from interf")

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
