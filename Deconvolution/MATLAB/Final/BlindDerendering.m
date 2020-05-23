function [s_est_out, brdf_est_out] = BlindDerendering(g, tau, tol, alpha, obs, doDisplay, mat, source, prop, maxCount)

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
    
    [s_est, brdf_est] = DerenderingOpt(g, s_est, brdf_est, T, nmove, nsource, tau, tol, alpha, wall_points, wall_points_ids, doDisplay, prop, maxCount);   
    
    reorder_vec = zeros(max(1, round(nmove/tau)),tau);
    for i=1:tau
        for j=1:max(1, round(nmove/tau))
            [~, I] = max(brdf_est(:, i, j));
            reorder_vec(j, i) = I;
        end
    end
    reorder_vec = mean(reorder_vec, 1);
    [~, idx] = sort(reorder_vec);
    
    brdf_est_out = zeros(size(brdf_est));
    s_est_out = zeros(size(s_est));
    
    for i=1:tau
        s_est_out(idx(i), :) = s_est(i, :);
        brdf_est_out(:, i, :) = brdf_est(:, idx(i), :);
    end

end

