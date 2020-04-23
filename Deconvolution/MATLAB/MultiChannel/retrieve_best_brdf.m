clear all
close all

iternum = 20;
pos = linspace(-100, -1, iternum);
err = zeros(1, iternum);
for i=1:iternum
    err(i) = iter(pos(i));
end
stem(pos, err)


%% FUNCTIONS

function out = iter(pos)
    num_lin = 50;
    num_ang = 360;
    margin = 5;
    mirror_BRDF = @(angle_diff, margin) (1*(abs(angle_diff) <= margin/2));

    sigma = 1;
    blurred_mirror_BRDF = @(angle_diff, sigma) (normpdf(-angle_diff, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 15, sigma)/normpdf(0, 0, sigma) + normpdf(-angle_diff, -35, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 35, sigma)/normpdf(0, 0, sigma));
    %blurred_mirror_BRDF = @(angle_diff, margin) (1*(abs(angle_diff+30) <= margin/2));

    obs_pos = [0, -3];
    obs_size = 4;
    obs_interval = [-45, 45];
    obs = build_obs(obs_pos, obs_size, obs_interval);

    source_pos = [0, pos];
    source_support_width = 10;
    source_support_size = 2;
    [source, interf_test] = build_source(source_pos, source_support_width, source_support_size, obs_size);

    [mirr_x, mirr_h, mirr_g, mirr_co] = rendering(obs, source, mirror_BRDF, num_lin, num_ang, margin);
    indices = mirr_g(1, :)>=1e-3;

    [x, h, g, co] = rendering(obs, source, blurred_mirror_BRDF, num_lin, num_ang, sigma);
    h_wall = compute_brdf_map_wall(xcorr(h, 'normalize'), obs, num_lin);
    %[x, h, g, co] = rendering(obs, source, mirror_BRDF, num_lin, num_ang, margin);
    [interf_brdf, interf_fbd, interf_g] = fbd(g, source_support_size, indices); 
    interf_brdf_true_wall = xcorr(h_wall);

%     figure;
%     subplot(1, 2, 1);
%     stem(interf_brdf);
%     subplot(1, 2, 2);
%     stem(h_wall);

    err = immse(interf_brdf, h_wall');
    out = err;

    % figure;
    % for i=1:obs(5)*obs(5)
    %     subplot(1, obs(5)*obs(5), i);
    %     stem(interf_g(:, i));
    % end
end

function [source, interf_test] = build_source(source_pos, source_support_width, source_support_size, obs_size)
    source = zeros(source_support_size, 3, obs_size);
    interf_test_prep = zeros(source_support_size, obs_size);
    if(source_support_size == 1)
        x_axis = zeros(1, 1);
    else
        x_axis = linspace(-source_support_width/2, source_support_width/2, source_support_size);
    end
    for i=1:obs_size
        for j=1:source_support_size  
            val = rand();
            source(j, :, i) = [source_pos(1)+x_axis(j), source_pos(2), val];
            interf_test_prep(j, i) = val;
        end
    end
    interf_test = xcorr(interf_test_prep);
end

function obs = build_obs(obs_pos, obs_size, obs_interval)
    obs = [obs_pos(1) obs_pos(2) obs_interval(1) obs_interval(2) obs_size];    
end

function intersection = ray_wall_intersection(ray_origin, ray_direction)
    v1 = ray_origin - [-1e20, 0];
    v2 = [1e20, 0];
    v3 = [-ray_direction(2), ray_direction(1)];
    dott = dot(v2, v3);
    if abs(dott) < 0.000001
        intersection = NaN;
        return;
    end
    t1 = (v2(1)*v1(2) - v2(2)*v1(1))/dott;
    t2 = dot(v1, v3)/dott;
    if t1 >= 0.0 && (t2 >= 0.0 && t2 <= 1.0)
        intersection = t1;
        return;
    end
    intersection = NaN;
    return;
    
end

% Convention : Wall center is 0, 0
function [x_axis, brdf_map, l_out, closest_obs] = rendering(obs, source, brdf, num_lin, num_ang, param)

    l_out = zeros(obs(5), num_lin);
    x_axis = zeros(obs(5), num_lin);
    brdf_map = zeros(1, num_ang);    
    closest_obs = 1;
    biggest_val = 0;    
    
    angles_brdf_map = linspace(-180, 180, num_ang);
    for i=1:length(angles_brdf_map)
        brdf_map(i) = brdf(angles_brdf_map(i), param);
    end
    
%     figure;
%     subplot(1, 2 + obs(5), 1);
%     stem(source(:,3));
%     subplot(1, 2 + obs(5), 2);
%     stem(angles_brdf_map, brdf_map);
    
    obs_angle = linspace(obs(3)+90, obs(4)+90, 2);
    wall_points = zeros(2, 2);    
    for j=1:2
        obs_direction = [cos(obs_angle(j)*pi/180), sin(obs_angle(j)*pi/180)];
        wall_points(j, :) = obs([1 2]) + obs_direction*ray_wall_intersection(obs([1 2]), obs_direction);
    end
    
    for i=1:obs(5)       
        x_axis(i, :) = linspace(wall_points(1, 1), wall_points(2, 1), num_lin);
        l_out_temp = zeros(1, num_lin);
        
        if source(1, 3, i) > biggest_val
            biggest_val = source(1, 3, i);
            closest_obs = i;
        end

        for j=1:round(num_lin)            
           
            normal = [x_axis(i, j), -1] - [x_axis(i, j), 0];
            dir_out = obs([1 2]) - [x_axis(i, j), 0];
            angle_out = atan2d(normal(1)*dir_out(2)-normal(2)*dir_out(1),normal(1)*dir_out(1)+normal(2)*dir_out(2));
            
            for k=1:size(source, 1)                
                dir_in = [source(k, 1, i), source(k, 2, i)] - [x_axis(i, j), 0];
                angle_in = atan2d(normal(1)*dir_in(2)-normal(2)*dir_in(1),normal(1)*dir_in(1)+normal(2)*dir_in(2));  
                l_out_temp(j) = l_out_temp(j) + source(k, 3, i)*brdf(angle_in + angle_out, param);                
            end           
            
        end
        
        l_out(i, :) = l_out_temp;
%         subplot(1, 2 + obs(5), 2 + obs(5)-i+1);
%         stem(x_axis(i, :), l_out(i, :));
    end   
    
end

function [s_interf_est, g_interf_est, d_interf] = fbd(d, tau, indices)
    
    % A few initializations
    indices_xcorr = xcorr(indices);
    indices_xcorr = indices_xcorr > 0.1;
    n = length(d(:, 1));
    T = length(d(1, :));
    % Prep FIBD : Make interferometric outputs
    d = d';
    d_interf = xcorr(d, 'normalize');
    
    % FIBD
    [s_interf_est, g_interf_est] = fibd(d_interf, tau, 1e-6, n, T, indices_xcorr);
end

function [s_interf_est, g_interf_est] = fibd(d_interf, tau, tol, n, T, indices)

    % Initialization
    s_interf_est = zeros(2*T-1, 1);
    s_interf_est(T) = 1;
    g_interf_est = rand(2*tau-1, n*n);
    for i = 1:n
        for j = 1:n
            index = get_col_num(i, j, n);
            if i == j 
                g_interf_est(:, index) = 0;
                g_interf_est(tau, index) = rand();
            end
        end
    end
    alpha = [1e5 0];    
    % Loop over decreasing alpha
    for i = 1:length(alpha)
        W1 = inf;
        W2 = inf;
        W1p = W1;
        W2p = W2;
        deltaW = inf; 
        count = 0;
        while deltaW > tol
            % 1. Update s
            % 1.1 Optimize W with s
            if(count > 15)
                break;
            end
            sa = optimvar('sa', 2*T-1);
            Wfun = @(sa) computeW(sa, g_interf_est, d_interf, n, tau, alpha(i), indices);
            Wexp = fcn2optimexpr(Wfun,sa);
            Wprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Wexp);
            Wprob.Constraints.cons1 = sa(T) == 1;            
            Wprob.Constraints.cons2 = sa(1:T-1) == flip(sa(T+1:end));
            Wprob.Constraints.cons3 = sa(:) >= 0; 
            W0.sa = s_interf_est;            
            [Wsol,Wfval,~,~] = solve(Wprob,W0,'Options', optimoptions(@fmincon, 'MaxFunctionEvaluations', 1e5));
            % 1.2 Assignements
            s_interf_est = Wsol.sa;
            W1p = W1;
            W1 = Wfval;
            % 2. Update g
            % 2.1 Optimize W with g            
            gij = optimvar('gij', 2*tau-1, n*n);
            Wfun = @(gij) computeW(s_interf_est, gij, d_interf, n, tau, alpha(i), indices);
            Wexp = fcn2optimexpr(Wfun,gij);
            Wprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Wexp); 
            W0.gij = g_interf_est;            
            [Wsol,Wfval,~,~] = solve(Wprob,W0,'Options', optimoptions(@fminunc, 'MaxFunctionEvaluations', 1e5));
            % 2.2 Assignements
            g_interf_est = Wsol.gij;  
            W2p = W2;
            W2 = Wfval;
            deltaW = max([W1p - W1 W2p - W2]);            
            disp(deltaW);
            count = count + 1;
        end
    end
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end

function out = computeW(s_interf_est, g_interf_est, d_interf, n, tau, alpha, indices)
    % 1.1.1 Compute V
    V = 0;
    for k=1:n
        for l=k:n
            term1 = abs(conv(s_interf_est, g_interf_est(:, get_col_num(k, l, n)), 'same'));
            term2 = d_interf(:, get_col_num(k, l, n));
            temp_diff = term2 - term1;
            temp_diff = temp_diff .* temp_diff;
            V = V + sum(temp_diff);
        end
    end
    % 1.1.2 Compute next term
    g_same = g_interf_est(:, get_col_num(1:n, 1:n, n));
    t = (-tau+1:tau-1)';
    Wterm = sum(sum(t.*t.*g_same.*g_same));
    W = V + alpha*Wterm;    
    out = W;
end

function out = fit_impulse(indices, source)
    ones_flag = 0;
    out = zeros(length(indices), 1);
    count = 1;
    for i=1:length(indices)
        if(indices(i)>=1)
            if(ones_flag == 0)
                out(i) = source(count);
                ones_flag = 1;
                count = count + 1;
            end      
        else
            ones_flag = 0;
        end
    end    
end

function out = compute_brdf_map_wall(brdf_map, obs, num_lin)
    % CREATE BRDF MAP
    obs_angle = linspace(obs(3)+90, obs(4)+90, 2);
    wall_points = zeros(2, 2);    
    for j=1:2
        obs_direction = [cos(obs_angle(j)*pi/180), sin(obs_angle(j)*pi/180)];
        wall_points(j, :) = obs([1 2]) + obs_direction*ray_wall_intersection(obs([1 2]), obs_direction);
    end
    angle_diff_vec = linspace(-180, 180, length(brdf_map));
    brdf_map_wall = zeros(1, 2*num_lin-1);
    x_axis = linspace(wall_points(1, 1), wall_points(2, 1), 2*num_lin-1);
    for j=1:length(brdf_map)
        ang = angle_diff_vec(j);
        wp = obs(2)*tand(ang);
        [~,idx] = min(abs(x_axis-wp));
        brdf_map_wall(idx) = brdf_map_wall(idx) + brdf_map(j);
    end    
    out = brdf_map_wall./max(brdf_map_wall);
end