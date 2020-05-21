function out = ComputeROpt(s_est, brdf_est_sparse, brdf_est_sparse_indices, brdf_est_dim1, brdf_est_dim2, g_flat, alpha, tau, nmove, nsource, wall_points, wall_points_ids, T)
    % 1.1.1 Compute V
    R = 0;
    brdf_est_opt = UnSparsify(brdf_est_sparse, brdf_est_sparse_indices, brdf_est_dim1, brdf_est_dim2);
    brdf_est_flat = BuildFlatFromOpt(brdf_est_opt, T, nmove, tau, wall_points, wall_points_ids);    
    for j=1:nmove
        g = g_flat(:, 1+(j-1)*nsource:j*nsource);
        brdf_est = brdf_est_flat(:, 1+(j-1)*tau:j*tau);
        temp = g - brdf_est*s_est;
        temp_diff = temp .* temp;
        R = R + sum(sum(temp_diff));
    end
    R = R + alpha*norm(brdf_est_opt, 1);
    out = R;
end
