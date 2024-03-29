clear all
close all

num_lin = 5000;
num_ang = 5000;
margin = 1;
mirror_BRDF = @(angle_diff, margin) (1*(abs(angle_diff) <= margin/2));

sigma = 1;
blurred_mirror_BRDF = @(angle_diff, sigma) (normpdf(-angle_diff, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 15, sigma)/normpdf(0, 0, sigma) + normpdf(-angle_diff, -35, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 35, sigma)/normpdf(0, 0, sigma));
%blurred_mirror_BRDF = @(angle_diff, margin) (1*(abs(angle_diff+30) <= margin/2));

obs_pos = [0, -5];
obs_size = 1;
obs_interval = [-45, 45];
obs = build_obs(obs_pos, obs_size, obs_interval);

source_pos = [0, -5];
source_support_width = 10;
source_support_size = 2;
[source, interf_test] = build_source(source_pos, source_support_width, source_support_size, obs_size);

[mirr_x, mirr_h, mirr_g] = rendering(obs, source, mirror_BRDF, num_lin, num_ang, margin, false, false);
[mirr_xc, mirr_hc, mirr_gc] = rendering(obs, source, mirror_BRDF, num_lin, num_ang, margin, true, false);

[x, h, g] = rendering(obs, source, blurred_mirror_BRDF, num_lin, num_ang, sigma, false, false);
[xc, hc, gc] = rendering(obs, source, blurred_mirror_BRDF, num_lin, num_ang, sigma, true, false);

err1 = immse(mirr_g, mirr_gc);
err2 = immse(g, gc);

figure;
subplot(1, 2, 1)
stem(mirr_g);
subplot(1, 2, 2)
stem(mirr_gc);

figure;
subplot(1, 2, 1)
stem(g);
subplot(1, 2, 2)
stem(gc);

% [mirr_x, mirr_h, mirr_g] = rendering(obs, source, mirror_BRDF, num_lin, num_ang, margin, false, false);
% mirr_g_xcorr = xcorr(mirr_g, 'normalized');
% [mirr_xc, mirr_hc, mirr_gc_xcorr] = rendering(obs, source, mirror_BRDF, 2*num_lin-1, num_ang, margin, true, true);
% 
% [x, h, g] = rendering(obs, source, blurred_mirror_BRDF, num_lin, num_ang, sigma, false, false);
% g_xcorr = xcorr(g', 'normalized');
% [xc, hc, gc_xcorr] = rendering(obs, source, blurred_mirror_BRDF, 2*num_lin-1, num_ang, sigma, true, true);

%err1 = immse(mirr_g_xcorr, mirr_gc_xcorr);
%err2 = immse(g_xcorr, gc_xcorr);

% figure;
% subplot(1, 2, 1)
% stem(mirr_g_xcorr);
% subplot(1, 2, 2)
% stem(mirr_gc_xcorr);
% 
% figure;
% subplot(1, 2, 1)
% stem(g_xcorr);
% subplot(1, 2, 2)
% stem(gc_xcorr);


%% Functions

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
function [x_axis, brdf_map, l_out] = rendering(obs, source, brdf, num_lin, num_ang, param, test_conv, test_xcorr)
    l_out = zeros(obs(5), num_lin);
    x_axis = zeros(obs(5), num_lin);
    %brdf_map = zeros(num_ang, num_lin);
    angles_brdf_map = linspace(-180, 180, num_ang).*ones(num_lin, num_ang);
    brdf_map = brdf(angles_brdf_map(1, :), param);
    %figure;
    %subplot(1, 2 + obs(5), 1);
    %stem(source(:,3));    
    
    obs_angle = linspace(obs(3)+90, obs(4)+90, 2);
    wall_points = zeros(2, 2);
    
    for j=1:2
        obs_direction = [cos(obs_angle(j)*pi/180), sin(obs_angle(j)*pi/180)];
        wall_points(j, :) = obs([1 2]) + obs_direction*ray_wall_intersection(obs([1 2]), obs_direction);
    end
    
    for i=1:obs(5)       
        
        if test_conv == false
            
            x_axis(i, :) = linspace(wall_points(1, 1), wall_points(2, 1), num_lin);
            l_out_temp = zeros(1, num_lin);
            for j=1:round(num_lin)           
                normal = [x_axis(i, j), -1] - [x_axis(i, j), 0];
                dir_out = obs([1 2]) - [x_axis(i, j), 0];
                angle_out = atan2d(normal(1)*dir_out(2)-normal(2)*dir_out(1),normal(1)*dir_out(1)+normal(2)*dir_out(2)); 
                %brdf_map(1:num_ang, j) = brdf(-(angles_brdf_map(1:num_ang) + angle_out), param);
                for k=1:size(source, 1)                
                    dir_in = [source(k, 1, i), source(k, 2, i)] - [x_axis(i, j), 0];
                    angle_in = atan2d(normal(1)*dir_in(2)-normal(2)*dir_in(1),normal(1)*dir_in(1)+normal(2)*dir_in(2));  
                    l_out_temp(j) = l_out_temp(j) + source(k, 3, i)*brdf(angle_in + angle_out, param);                
                end 
            end
            
        else        
            
            x_axis(i, :) = linspace(wall_points(1, 1), wall_points(2, 1), num_lin);
            l_out_temp = zeros(1, num_lin);            
            normal = [x_axis(i, 1:num_lin); -1*ones(1, num_lin)] - [x_axis(i, 1:num_lin); zeros(1, num_lin)];
            dir_out = [obs(1)*ones(1, num_lin); obs(2)*ones(1, num_lin)] - [x_axis(i, 1:num_lin); zeros(1, num_lin)];
            angle_out = atan2d(normal(1, :).*dir_out(2, :)-normal(2, :).*dir_out(1, :), normal(1, :).*dir_out(1, :)+normal(2, :).*dir_out(2, :));                 
            lumin_map = zeros(num_lin, num_ang);
            
            for k=1:size(source, 1)
                dir_in = [source(k, 1, i)*ones(1, num_lin); source(k, 2, i)*ones(1, num_lin)] - [x_axis(i, 1:num_lin); zeros(1, num_lin)];
                angle_in = atan2d(normal(1, :).*dir_in(2, :)-normal(2, :).*dir_in(1, :) ,normal(1, :).*dir_in(1, :)+normal(2, :).*dir_in(2, :));
                [~,idx] = min(abs(angles_brdf_map-angle_in'), [], 2);
                indexes = sub2ind(size(lumin_map),(1:num_lin)',idx);
                lumin_map(indexes) = lumin_map(indexes) + source(k, 3, i);
            end
            
            if(test_xcorr)
                prop = 0.1;
                temp1 = xcorr(lumin_map);
                temp2 = xcorr(brdf_map);
                temp = conv(temp1, temp2);
                temp = temp(round(prop*length(temp)):round(end-prop*length(temp)))./max(temp);
                %temp = temp./max(temp);
            else
                temp = zeros(num_lin, num_ang);
                for j=1:num_lin
                    temp(j, :) = conv(lumin_map(j,:), brdf_map, 'same'); % Make full BRDF map to avoid 'same'
                end                
            end      
            
            angles_temp = linspace(angles_brdf_map(1, 1), angles_brdf_map(1, end), num_ang).*ones(num_lin, num_ang);
            [~,idx] = min(abs(angles_temp + angle_out'), [], 2);
            indexes = sub2ind(size(angles_temp),(1:num_lin)',idx);
            l_out_temp(end-(1:num_lin)+1) = temp(indexes);
            
        end
        l_out(i, :) = l_out_temp;   
        
%         subplot(1, 2 + obs(5), 2);
%         imagesc(brdf_map);
%         
%         subplot(1, 2 + obs(5), 2 + obs(5)-i+1);        
        if test_conv == false
            %stem(x_axis(i, :), l_out(i, :));
        else
            l_out(i, :) = flip(l_out(i, :));
            %stem(x_axis(i, :), l_out(i, :));            
        end
        
    end   
    
end