function [errx, errh] = DerenderingWrapper(num_lin, obs_size_move, obs_size_source, source_support_size, prop, tol, alpha)
    sigma = 1;
    blurred_mirror_BRDF = @(angle_diff, sigma) (normpdf(-angle_diff, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 15, sigma)/normpdf(0, 0, sigma) + normpdf(-angle_diff, -35, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 35, sigma)/normpdf(0, 0, sigma));

    obs_pos = [0, -5];
    obs_interval = [-45, 45];
    obs = build_obs(obs_pos, obs_size_move, obs_size_source, obs_interval);

    source_pos = [0, -5];
    source_support_width = 10;
    [source, interf_test] = build_source(source_pos, source_support_width, source_support_size, obs_size_source, false);

    [x_axis, mat] = CreateRenderingMatrixFromBRDFMoveCam(obs, source, blurred_mirror_BRDF, num_lin, sigma, prop);
    g = zeros(num_lin, obs_size_move, obs_size_source);

    for i=1:obs_size_source
        for j=1:obs_size_move
            A = mat(:, :, j);
            x = source(:, 3, i);
            y = A*x;
            g(:, j, i) = y;       
        end
    end

    disp("Rendering done.");

    disp("Starting de-rendering");
    [x_est, h_est, h_est_flat] = blind_derendering(g, source_support_size, tol, alpha, false);

    errx = immse(reshape(source(:, 3, :),[source_support_size, obs_size_source]), x_est);
    errh = immse(h_est, mat);

end

function [x_est, h_est, h_est_flat] = blind_derendering(d_full, tau, tol, alpha, useParallel)
    T = length(d_full(:, 1, 1));
    nmove = length(d_full(1, :, 1));
    nsource = length(d_full(1, 1, :));
    g_est = zeros(T, tau, nmove);
    for i=1:tau
        for j=1:nmove
            mm = mean(d_full, 3);
            g_est(:, i, j) = mm(:, j);
        end
    end    
    s_est_full = rand(tau, nsource);

    
    [s_est_full, g_est, g_est_flat] = DerenderingStandardMatrixMoveCamGeneticComp(d_full, s_est_full, g_est, T, nmove, nsource, tau, tol, alpha, useParallel);   
    
    h_est = g_est;    
    h_est_flat = g_est_flat;
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