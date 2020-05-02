clear all
close all

num_lin = 4;
num_ang = 25;
margin = 2;
mirror_BRDF = @(angle_diff, margin) (1*(abs(angle_diff) <= margin/2));
sigma = 1;
blurred_mirror_BRDF = @(angle_diff, sigma) (normpdf(-angle_diff, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 15, sigma)/normpdf(0, 0, sigma) + normpdf(-angle_diff, -35, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 35, sigma)/normpdf(0, 0, sigma));

obs_pos = [0, -5];
obs_size_move = 40;
obs_size_source = 40;
obs_interval = [-45, 45];
obs = build_obs(obs_pos, obs_size_move, obs_size_source, obs_interval);

source_pos = [0, -5];
source_support_width = 10;
source_support_size = 4;
[source, interf_test] = build_source(source_pos, source_support_width, source_support_size, obs_size_source, false);
[empty_source, ~] = build_source(source_pos, source_support_width, source_support_size, 1, true);

[x_axis, mat] = CreateRenderingMatrixFromBRDFMoveCam(obs, source, blurred_mirror_BRDF, num_lin, sigma, 0.5);
g = zeros(num_lin, obs_size_move, obs_size_source);

for i=1:obs_size_source
    %figure;
    for j=1:obs_size_move
        A = mat(:, :, j);
        x = source(:, 3, i);
        y = A*x;
        g(:, j, i) = y;
        %subplot(1, obs_size_move, j);
        %stem(flip(x_axis(:, obs_size_move-j+1)), g(:, j, i));        
    end
end
g_test = reshape(g(:, :, 1), [num_lin, obs_size_move]);

% obs(5) = 1;
% [x, h, g_verif] = FullRendering(obs, source, blurred_mirror_BRDF, num_lin, num_ang, sigma);
% figure;
% stem(g_verif);

[x_est, h_est] = blind_derendering(g, source_support_size, source, mat);
figure;
for i=1:obs_size_move
    subplot(2, obs_size_move, i);
    imagesc(h_est(:, :, i));
    subplot(2, obs_size_move, i+obs_size_move);
    imagesc(mat(:, :, i));
end


function [x_est, h_est] = blind_derendering(d_full, tau, source, mat)
    T = length(d_full(:, 1, 1));
    nmove = length(d_full(1, :, 1));
    nsource = length(d_full(1, 1, :));
    mat_flat = reshape(mat, [T, tau*nmove]);
    %d_interf = xcorr(d, 'normalize');
    d_interf = zeros(2*T-1, nmove*nmove, nsource);
    for i=1:nsource
        d_interf(:, :, i) = xcorr(d_full(:, :, i));
    end
    g_est = zeros(T, tau, nmove);
    for i=1:tau
        for j=1:nmove
            mm = mean(d_full, 3);
            g_est(:, i, j) = mm(:, j);
        end
    end
    
    s_est_full = rand(tau, nsource);
    %s_est_full = reshape(source(:, 3, :), [tau, nsource]);
    %s_est(round(tau/2)) = 1;
    s_interf_est = rand(2*T-1, tau);
    %s_interf_est = mat_interf;
    
    %[g_cov, s_est_full] = EstimateSourceMoveCam(d_full, s_est_full, g_est, T, nmove, nsource, tau, 1e-8);   
    [s_est_full, g_est, g_est_flat] = DerenderingStandardMatrixMoveCamGenetic(d_full, s_est_full, mat, T, nmove, nsource, tau, 1e-3);
    %s_est = s_est./max(s_est);
    
    figure;
    subplot(1, 2, 1);
    imagesc(mat_flat);
    subplot(1, 2, 2);
    imagesc(g_est_flat);
    
    figure;
    subplot(1, 2, 1);
    imagesc(s_est_full);
    subplot(1, 2, 2);
    imagesc(reshape(source(:, 3, :), [tau, nsource]));
    
    h_est = g_est;    
    x_est = s_est_full;
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

function obs = build_obs(obs_pos, obs_size_move, obs_size_source, obs_interval)
    obs = [obs_pos(1) obs_pos(2) obs_interval(1) obs_interval(2) obs_size_move obs_size_source];    
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end