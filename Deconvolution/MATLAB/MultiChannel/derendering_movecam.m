clear all
close all

num_lin = 25;
num_ang = 25;
margin = 2;
mirror_BRDF = @(angle_diff, margin) (1*(abs(angle_diff) <= margin/2));
sigma = 1;
blurred_mirror_BRDF = @(angle_diff, sigma) (normpdf(-angle_diff, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 15, sigma)/normpdf(0, 0, sigma) + normpdf(-angle_diff, -35, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 35, sigma)/normpdf(0, 0, sigma));

obs_pos = [0, -5];
obs_size = 4;
obs_interval = [-45, 45];
obs = build_obs(obs_pos, obs_size, obs_interval);

source_pos = [0, -5];
source_support_width = 10;
source_support_size = 4;
[source, interf_test] = build_source(source_pos, source_support_width, source_support_size, 1, false);
[empty_source, ~] = build_source(source_pos, source_support_width, source_support_size, 1, true);

[x_axis, rendering_matrices] = CreateRenderingMatrixFromBRDFMoveCam(obs, source, blurred_mirror_BRDF, num_lin, sigma, 0.5);
g = zeros(num_lin, obs_size);
figure;
for i=1:obs_size
    g(:, i) = rendering_matrices(:, :, i)*source(:, 3, 1);
    subplot(1, obs_size, i);
    stem(flip(x_axis(:, obs_size-i+1)), g(:, i));
end
rendering_matrices_sum = reshape(sum(rendering_matrices, 2)/source_support_size, num_lin, obs_size);


% obs(5) = 1;
% [x, h, g_verif] = FullRendering(obs, source, blurred_mirror_BRDF, num_lin, num_ang, sigma);
% figure;
% stem(g_verif);

[x_est, h_est, x_interf_est, h_interf_est] = blind_derendering(g, source_support_size, obs, empty_source, num_ang, num_lin, mirror_BRDF, margin, rendering_matrices_sum, xcorr(rendering_matrices_sum), xcorr(g), source, rendering_matrices);
figure;
for i=1:obs_size
    subplot(2, obs_size, i);
    imagesc(h_est(:, :, i));
    subplot(2, obs_size, i+obs_size);
    imagesc(rendering_matrices(:, :, i));
end

function [x_est, h_est, x_interf_est, h_interf_est] = blind_derendering(d, tau, obs, empty_source, num_ang, num_lin, mirror_BRDF, margin, mat_interf1, mat_interf2, mat_interf3, source, rendering_matrices)
    n = length(d(1, :));
    T = length(d(:, 1));
    %d_interf = xcorr(d, 'normalize');
    d_interf = xcorr(d);
    g_est = rand(T, tau, n);
    s_est = ones(tau, 1);
    %s_est = source(:, 3, 1);
    %s_est(round(tau/2)) = 1;
    s_interf_est = rand(2*T-1, tau);
    %s_interf_est = mat_interf;
    
    g_interf_est = RecovXcorrMovecam(d, d_interf, n, T, tau);
    g_est_summed = PhaseRetrieval(g_interf_est, T, n, 1);
    [s_est, g_est] = DerenderingStandardMatrixMoveCam(d, s_est, g_est_summed, g_est, T, n, tau, obs, empty_source, num_ang, 1e-8);
    %s_est = s_est./max(s_est);
    h_est = g_est;
    h_interf_est = s_interf_est;
    x_est = s_est;
    x_interf_est = g_interf_est;
end

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
%         if ~isempty
%             source(:, 3, i) = 0.9*source(:, 3, i);
%             ind = randi(source_support_size);
%             source(ind, 3, i) = 1;
%         end
    end
    interf_test = xcorr(interf_test_prep);
end

function obs = build_obs(obs_pos, obs_size, obs_interval)
    obs = [obs_pos(1) obs_pos(2) obs_interval(1) obs_interval(2) obs_size];    
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end


