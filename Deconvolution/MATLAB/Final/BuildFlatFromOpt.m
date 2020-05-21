function out = BuildFlatFromOpt(brdf_est_opt, T, nmove, tau, wall_points, wall_points_ids)

    out = zeros(T, tau*nmove);
    init = round(nmove/2);
    
    for j=init:nmove
        wp = wall_points(j, :);
        wp_ids = wall_points_ids(j, :);
        if j==init
            out(:, 1+(j-1)*tau:j*tau) = brdf_est_opt(1:T, :);
            g_prev = brdf_est_opt(1:T, :);
            curr_ind = T+1;
        else            
            num = sum(wp > prev_wp(end));            
            g_est = zeros(T, tau);
            g_est(1:end-num, :) = g_prev(wp_ids(1:end-num), :);
            g_est(end-num+1:end, :) = brdf_est_opt(curr_ind:curr_ind+num-1, :);
            out(:, 1+(j-1)*tau:j*tau) = g_est;
            g_prev = g_est;
            curr_ind = curr_ind+num;
        end
        prev_wp = wp;
    end
    
    g_prev = brdf_est_opt(1:T, :);
    prev_wp = wall_points(init, :);
    for j=init-1:-1:1
        wp = wall_points(j, :);  
        wp_ids = wall_points_ids(j, :);
        num = sum(wp < prev_wp(1));        
        g_est = zeros(T, tau);
        g_est(num+1:end, :) = g_prev(wp_ids(num+1:end), :);
        g_est(1:num, :) = brdf_est_opt(curr_ind:curr_ind+num-1, :);
        out(:, 1+(j-1)*tau:j*tau) = g_est;
        g_prev = g_est;
        curr_ind = curr_ind+num;
        prev_wp = wp;
    end
end

