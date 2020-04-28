clear all
%close all

num_lin = 2000;
num_ang = 2000;
margin = 5;
mirror_BRDF = @(angle_diff, margin) (1*(abs(angle_diff) <= margin/2));

sigma = 1;
blurred_mirror_BRDF = @(angle_diff, sigma) (normpdf(-angle_diff, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 15, sigma)/normpdf(0, 0, sigma));
%blurred_mirror_BRDF = @(angle_diff, margin) (1*(abs(angle_diff+30) <= margin/2));

obs_pos = [0, -5];
obs_size = 4;
gap = 1/obs_size;
obs_interval = [-45, 45];
obs = build_obs(obs_pos, obs_size, obs_interval);

source_pos = [0, -5];
source_support_width = 10;
source_support_size = 2;
source_function = @(x_axis) (0.25*(x_axis < 0) + (1*x_axis >= 0));
%source_function = @(x_axis) (1);
source = build_source(source_pos, source_support_width, source_support_size, source_function);

[mirr_x, mirr_h, mirr_g, mirr_co] = rendering(obs, source, mirror_BRDF, num_lin, num_ang, margin, gap);
[x, h, g, co] = rendering(obs, source, blurred_mirror_BRDF, num_lin, num_ang, sigma, gap);
% [h_est, x_est] = fbd(g, 2, co); 
% figure;
% subplot(1, 2, 1);
% stem(h_est);
% subplot(1, 2, 2);
% stem(xcorr(h));
% figure;
% for i=1:obs_size
%     subplot(2, obs_size, i);
%     stem(x_est(:,i));
% end
% for i=1:obs_size
%     subplot(2, obs_size, obs_size + i);
%     stem(mirr_g(i,:));
% end

function source = build_source(source_pos, source_support_width, source_support_size, source_function)
    source = zeros(source_support_size, 3);
    if(source_support_size == 1)
        x_axis = zeros(1, 1);
    else
        x_axis = linspace(-source_support_width/2, source_support_width/2, source_support_size);
    end    
    for i=1:source_support_size
        source(i, :) = [source_pos(1)+x_axis(i), source_pos(2), source_function(x_axis(i))];
    end
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
function [x_axis, brdf_map, l_out, closest_obs] = rendering(obs, source, brdf, num_lin, num_ang, param, gap)

    l_out = zeros(obs(5), round(num_lin/obs(5)));
    x_axis = zeros(obs(5), round(num_lin/obs(5)));
    brdf_map = zeros(1, num_ang);    
    closest_obs = 1;
    smallest_dist = inf;    
    
    angles_brdf_map = linspace(-180, 180, num_ang);
    for i=1:length(angles_brdf_map)
        brdf_map(i) = brdf(angles_brdf_map(i), param);
    end
    
    figure;
    subplot(1, 2 + obs(5), 1);
    stem(source(:,3));
    subplot(1, 2 + obs(5), 2);
    stem(angles_brdf_map, brdf_map);
    
    obs_angle = linspace(obs(3)+90, obs(4)+90, 2);
    wall_points = zeros(2, 2);
    for j=1:2
        obs_direction = [cos(obs_angle(j)*pi/180), sin(obs_angle(j)*pi/180)];
        wall_points(j, :) = obs([1 2]) + obs_direction*ray_wall_intersection(obs([1 2]), obs_direction);
    end
    iter = wall_points(2, 1) - wall_points(1, 1) - gap;
    current_x_left = wall_points(1, 1);
    current_x_right = current_x_left + iter;
    
    for i=1:obs(5)       
        x_axis(i, :) = linspace(current_x_left, current_x_right, round(num_lin/obs(5)));
        l_out_temp = zeros(1, round(num_lin/obs(5)));        

        for j=1:round(num_lin/obs(5)) 
            
            dist_temp = norm([obs(1) obs(2)]-[x_axis(i, j) 0]);
            if dist_temp < smallest_dist
                closest_obs = i;
            end
            
            normal = [x_axis(i, j), -1] - [x_axis(i, j), 0];
            dir_out = obs([1 2]) - [x_axis(i, j), 0];
            angle_out = atan2d(normal(1)*dir_out(2)-normal(2)*dir_out(1),normal(1)*dir_out(1)+normal(2)*dir_out(2));
            
            for k=1:size(source, 1)                
                dir_in = [source(k, 1), source(k, 2)] - [x_axis(i, j), 0];
                angle_in = atan2d(normal(1)*dir_in(2)-normal(2)*dir_in(1),normal(1)*dir_in(1)+normal(2)*dir_in(2));  
                l_out_temp(j) = l_out_temp(j) + source(k, 3)*brdf(angle_in + angle_out, param);                
            end           
            
        end 
        
        l_out(i, :) = l_out_temp;      
        current_x_left = current_x_left + gap;
        current_x_right = current_x_left + iter;
        
        subplot(1, 2 + obs(5), 2 + obs(5)-i+1);
        stem(x_axis(i, :), l_out(i, :));
    end   
    
end

function [s_est, g_est, g_interf_est] = fbd(d, tau, closest_obs)
    
    % A few initializations
    n = length(d(:, 1));
    T = length(d(1, :));
    % Prep FIBD : Make interferometric outputs
    d = d';
    d_interf = xcorr(d, 'normalize');
    
    % FIBD
    [s_interf_est, g_interf_est] = fibd(d_interf, tau, 1e-5, n, T);
    g_est =  fpr(g_interf_est, tau, n, 10); 
    s_est = s_interf_est;
end

function [s_interf_est, g_interf_est] = fibd(d_interf, tau, tol, n, T)

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
        while deltaW > tol
            % 1. Update s
            % 1.1 Optimize W with s
            sa = optimvar('sa', 2*T-1);
            Wfun = @(sa) computeW(sa, g_interf_est, d_interf, n, tau, alpha(i));
            Wexp = fcn2optimexpr(Wfun,sa);
            Wprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Wexp);
            Wprob.Constraints.cons1 = sa(T) == 1;            
            Wprob.Constraints.cons2 = sa(1:T-1) == flip(sa(T+1:end));
            %Wprob.Constraints.cons3 = sa(1:end) >= 0;
            W0.sa = s_interf_est;            
            [Wsol,Wfval,~,~] = solve(Wprob,W0,'Options', optimoptions(@fmincon,'Display','iter', 'MaxFunctionEvaluations', 1e5));
            % 1.2 Assignements
            s_interf_est = Wsol.sa;
            W1p = W1;
            W1 = Wfval;
            % 2. Update g
            % 2.1 Optimize W with g
            gij = optimvar('gij', 2*tau-1, n*n);
            Wfun = @(gij) computeW(s_interf_est, gij, d_interf, n, tau, alpha(i));
            Wexp = fcn2optimexpr(Wfun,gij);
            Wprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Wexp); 
            %Wprob.Constraints.cons1 = gij(:) >= 0.1;
            %Wprob.Constraints.cons2 = gij(:) <= tau;
            W0.gij = g_interf_est;            
            [Wsol,Wfval,~,~] = solve(Wprob,W0,'Options', optimoptions(@fminunc,'Display','iter', 'MaxFunctionEvaluations', 1e5));
            % 2.2 Assignements
            g_interf_est = Wsol.gij;
            W2p = W2;
            W2 = Wfval;
            deltaW = max([W1p - W1 W2p - W2]);            
            disp(deltaW);
        end
    end
end

function g_est = fpr(g_interf_est, tau, n, f)
    
    beta = [1e5 0];
    g_est = rand(tau, n);
    g_est(1, f) = rand();
    g_est(2:end, f) = 0;
    
    for i = 1:length(beta)       
        % 1. Update g
        % 1.1 Optimize W with s
        gi = optimvar('gi', tau, n);
        Yfun = @(gi) computeY(n, g_interf_est, gi, tau, f, beta(i));
        Yexp = fcn2optimexpr(Yfun,gi);
        Yprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Yexp);
        Y0.gi = g_est;            
        [Ysol,Yfval,~,~] = solve(Yprob,Y0,'Options', optimoptions(@fminunc,'Display','iter', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e6));
        g_est = Ysol.gi;
    end
    gi = optimvar('gi', tau, n);
    Xfun = @(gi) computeX(n, g_interf_est, gi);
    Xexp = fcn2optimexpr(Xfun,gi);
    Xprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Xexp);
    %Xprob.Constraints.cons1 = gi(:) >= 0;
    X0.gi = g_est;            
    [Xsol,Xfval,~,~] = solve(Xprob,X0,'Options', optimoptions(@fminunc,'Display','iter', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e6));
    g_est = Xsol.gi;
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end

function out = computeW(s_interf_est, g_interf_est, d_interf, n, tau, alpha)
    % 1.1.1 Compute V
    V = 0;
    for k=1:n
        for l=k:n
            temp_diff = d_interf(:, get_col_num(k, l, n)) - conv(s_interf_est, g_interf_est(:, get_col_num(k, l, n)), 'same');
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

function out = computeX(n, g_interf_est, g_est)
    X = 0;
    for k=1:n
        for l=k:n
           temp_diff = g_interf_est(:, get_col_num(k, l, n)) - xcorr(g_est(:,k), g_est(:,l));
           temp_diff = temp_diff .* temp_diff;
           X = X + sum(temp_diff);
        end
    end
    out = X;
end

function out = computeY(n, g_interf_est, g_est, tau, f, beta)
    Y1 = 0;
    for k=1:n
       temp_diff = g_interf_est(:, get_col_num(k, f, n)) - xcorr(g_est(:,k), g_est(:,f));
       temp_diff = temp_diff .* temp_diff;
       Y1 = Y1 + sum(temp_diff);        
    end
    t =(1:tau).*(1:tau);
    gf = g_est(:,f).*g_est(:,f);
    Y2 = beta.*sum(gf.*t');
    out = Y1 + Y2;
end
