function out = BuildOptFromFlat(brdf_est_flat, T, nmove, tau, wall_points)    
    out = zeros(T*nmove, tau);
    init = round(nmove/2);
    for j=init:nmove
        wp = wall_points(j, :);
        brdf_est = brdf_est_flat(:, 1+(j-1)*tau:j*tau);
        if j==init
            out(1:T, :) = brdf_est;
            curr_ind = T+1;
        else
            num = sum(wp > prev_wp(end));
            out(curr_ind:curr_ind+num-1, :) = brdf_est(end-num+1:end, :);
            curr_ind = curr_ind+num;
        end
        prev_wp = wp;
    end
    prev_wp = wall_points(init, :);
    for j=init-1:-1:1
        wp = wall_points(j, :);
        brdf_est = brdf_est_flat(:, 1+(j-1)*tau:j*tau);
        num = sum(wp < prev_wp(1));
        out(curr_ind:curr_ind+num-1, :) = brdf_est(1:num, :);
        curr_ind = curr_ind+num;
        prev_wp = wp;
    end
    if nmove > 1
            out = out(1:curr_ind, :);
    end
end
