clear all
close all

num_lin = input("Number of wall sample points : ");
sigma = 1;
blurred_mirror_BRDF = @(angle_diff, sigma) (normpdf(-angle_diff, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 15, sigma)/normpdf(0, 0, sigma) + normpdf(-angle_diff, -35, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 35, sigma)/normpdf(0, 0, sigma));

obs_pos = [0, -5];
obs_size_move = input("Number of camera moves : ");
obs_size_source = input("Number of sources : ");
obs_interval = [-45, 45];
obs = build_obs(obs_pos, obs_size_move, obs_size_source, obs_interval);

source_pos = [0, -5];
source_support_width = 10;
source_support_size = input("Number of points by source : ");
[source, interf_test] = build_source(source_pos, source_support_width, source_support_size, obs_size_source, false);
x_true = reshape(source(:, 3, :),[source_support_size, obs_size_source]);

prop = input("Move proportion (default : 0.5) : ");
[x_axis, mat] = CreateRenderingMatrixFromBRDFMoveCam(obs, source, blurred_mirror_BRDF, num_lin, sigma, prop);
g = zeros(num_lin, obs_size_source, obs_size_move);

x = reshape(source(:, 3, :),[source_support_size, obs_size_source]);
for j=1:obs_size_move
    A = mat(:, :, j);    
    y = A*x;
    g(:, :, j) = y;       
end

disp("Rendering done.");

tol = input("Covergence tolerance : ");
alpha = input("L1 regularization coefficient (default : 1) : ");
beta = input("Autocorr regularization coefficient (default : 0.01) : ");
useParallel = input("Use parallel optimization (default : false) : ");

disp("Starting de-rendering");
[x_est, h_est, h_est_flat] = blind_derendering(g, source_support_size, tol, alpha, beta, useParallel, 'iter');
figure;
for i=1:obs_size_move
    subplot(2, obs_size_move, i);
    imagesc(h_est(:, :, i));
    subplot(2, obs_size_move, i+obs_size_move);
    imagesc(mat(:, :, i));
end

disp("De-rendering done. Writing output :");
save("derendering.mat")

function [x_est, h_est, h_est_flat] = blind_derendering(d_full, tau, tol, alpha, beta, useParallel, doDisplay)
    T = length(d_full(:, 1, 1));
    nmove = length(d_full(1, 1, :));
    nsource = length(d_full(1, :, 1));
    g_est = zeros(T, tau, nmove);
    for i=1:tau
        for j=1:nmove
            mm = mean(d_full, 2);
            g_est(:, i, j) = mm(:, j);
        end
    end    
    s_est_full = rand(tau, nsource);

    
    [s_est_full, g_est, g_est_flat] = DerenderingStandardMatrixMoveCamCascade(d_full, s_est_full, g_est, T, nmove, nsource, tau, tol, alpha, beta, useParallel, doDisplay);   
    
    reorder_vec = zeros(1,tau);
    for i=1:tau
        reorder_vec(i) = mean(g_est(:, i, 1), 1);
    end
    [~, idx] = sort(reorder_vec);
    
    h_est = zeros(size(g_est));
    h_est_flat = zeros(size(g_est_flat));
    x_est = zeros(size(s_est_full));
    
    for i=1:tau
        x_est(idx(i), :) = s_est_full(i, :);
        h_est_flat(:, idx(i):tau:end) = g_est_flat(:, i:tau:end);
        h_est(:, idx(i), :) = g_est(:, i, :);
    end

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