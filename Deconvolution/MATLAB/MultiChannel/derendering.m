clear all
close all

num_lin = 50;
num_ang = 50;
margin = 2;
mirror_BRDF = @(angle_diff, margin) (1*(abs(angle_diff) <= margin/2));
sigma = 1;
blurred_mirror_BRDF = @(angle_diff, sigma) (normpdf(-angle_diff, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 15, sigma)/normpdf(0, 0, sigma) + normpdf(-angle_diff, -35, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 35, sigma)/normpdf(0, 0, sigma));

obs_pos = [0, -5];
obs_size = 4;
obs_interval = [-45, 45];
obs = build_obs(obs_pos, obs_size, obs_interval);

source_pos = [0, -50];
source_support_width = 10;
source_support_size = 4;
[source, interf_test] = build_source(source_pos, source_support_width, source_support_size, obs_size, false);
[empty_source, ~] = build_source(source_pos, source_support_width, source_support_size, obs_size, true);

[x, h, g] = FullRendering(obs, source, blurred_mirror_BRDF, num_lin, num_ang, sigma);
mat_interf = (CreateRenderingMatrixFromBRDFInterf(obs, source, blurred_mirror_BRDF, num_lin, sigma))';
[x_est, h_est, x_interf_est, h_interf_est] = blind_derendering(g, source_support_size, obs, empty_source, num_ang, num_lin, mirror_BRDF, margin, mat_interf);

function [x_est, h_est, x_interf_est, h_interf_est] = blind_derendering(d, tau, obs, empty_source, num_ang, num_lin, mirror_BRDF, margin, mat_interf)
    n = length(d(:, 1));
    T = length(d(1, :));
    d = d';
    d_interf = xcorr(d, 'normalize');
    
    g_est = rand(tau, n);
    s_est = zeros(T, 1);
    s_est(round(T/2)) = 1;
    s_interf_est = rand(2*T-1, tau);
    s_interf_est = mat_interf;
    
    g_interf_est = rand(2*tau-1, n*n);
    
    %[~, g_interf_est] = FIBD(d_interf, T, n, tau, 5e-3, obs, empty_source, mirror_BRDF, num_lin, margin);
    [s_interf_est, g_interf_est] = DerenderingInterfMatrix(d_interf, s_interf_est, g_interf_est, T, n, tau, obs, empty_source, num_ang, 1e-8);
    g_est = PhaseRetrieval(g_interf_est, tau, n, false);
    s_est = PhaseRetrievalAuto(s_interf_est, T, tau);
    [s_est, g_est] = DerenderingPhaseRetrievalMatrix(d_interf, s_est, g_est, T, n, tau, obs, empty_source, num_ang, 1e-3);
    [s_est, g_est] = DerenderingStandardMatrix(d, s_est, g_est, T, n, tau, obs, empty_source, num_ang, 1e-8);
    for i=1:n
        g_est(:, i) = g_est(:, i)./max(g_est(:, i));
    end
    h_est = s_est;
    h_interf_est = s_interf_est;
    x_est = g_est;
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

