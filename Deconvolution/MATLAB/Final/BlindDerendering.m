function [s_est_out, brdf_est_out] = BlindDerendering(g, tau, tol, alpha, obs, doDisplay, mat, source, prop)

    T = length(g(:, 1, 1));
    nmove = length(g(1, 1, :));
    nsource = length(g(1, :, 1));
    brdf_est = zeros(T, tau, nmove);
    
    for i=1:tau
        for j=1:nmove
            mm = mean(g, 2);
            brdf_est(:, i, j) = mm(:, j);
        end
    end
    s_est = rand(tau, nsource);
    
    wall_points = ComputeWallPoints(obs, T, nmove);
    wall_points_ids = ComputeWallPointsIds(wall_points);
    
    [s_est, brdf_est] = DerenderingOpt(g, s_est, brdf_est, T, nmove, nsource, tau, tol, alpha, wall_points, wall_points_ids, doDisplay, prop);   
    
    reorder_vec = zeros(1,tau);
    for i=1:tau
        reorder_vec(i) = mean(brdf_est(:, i, 1), 1);
    end
    [~, idx] = sort(reorder_vec);
    
    brdf_est_out = zeros(size(brdf_est));
    s_est_out = zeros(size(s_est));
    
    for i=1:tau
        s_est_out(idx(i), :) = s_est(i, :);
        brdf_est_out(:, idx(i), :) = brdf_est(:, i, :);
    end

end

