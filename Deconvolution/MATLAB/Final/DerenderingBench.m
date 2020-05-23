function [s_mse, brdf_mse] = DerenderingBench(num_lin, angle, sigma, obs_size_move, obs_size_source, source_support_size, alpha, prop, max_count, tol)

    blurred_mirror_BRDF = @(angle_in, angle_out, sigma) (normpdf(angle_in+angle_out, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(angle_in+angle_out, 15, sigma)/normpdf(0, 0, sigma) + normpdf(angle_in+angle_out, -35, sigma)/normpdf(0, 0, sigma)  +  normpdf(angle_in+angle_out, 35, sigma)/normpdf(0, 0, sigma));

    obs_pos = [0, -5];
    obs_interval = [45, -45]; % Convention : Move from left to right
    obs_width = angle; % Must be smaller than the interval
    focus_dist = 2.5;
    obs = BuildObs(obs_pos, obs_size_move, obs_size_source, obs_interval, obs_width, focus_dist);

    source_pos = [0, -1];
    source_support_width = 5;
    source = BuildSource(source_pos, source_support_width, source_support_size, obs_size_source, false);
    s = reshape(source(:, 3, :),[source_support_size, obs_size_source]);

    [~, brdf_mat] = CreateRenderingMatrixFromBRDFFinal(obs, source, blurred_mirror_BRDF, num_lin, sigma);
    g = RenderingMatrixProduct(s, brdf_mat, obs_size_move, obs_size_source, num_lin);   

    [s_est, brdf_est] = BlindDerendering(g, source_support_size, tol, alpha, obs, 'off', brdf_mat, s, prop, max_count);
    brdf_mat_flat = Flatten(brdf_mat);
    brdf_est_flat = Flatten(brdf_est);
    
    brdf_mse = immse(brdf_mat_flat, brdf_est_flat);
    s_mse = immse(s, s_est);

end

