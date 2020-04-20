clear all
close all

step_lin = 1e-03;
step_ang = 1e-03;
margin = 1e-01;
mirror_BRDF = @(angle_diff, margin) (1*(abs(angle_diff) <= margin/2));

sigma = 1;
blurred_mirror_BRDF = @(angle_diff, sigma) (normpdf(-angle_diff, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(-angle_diff, 15, sigma)/normpdf(0, 0, sigma));
%blurred_mirror_BRDF = @(angle_diff, margin) (1*(abs(angle_diff+30) <= margin/2));


obs_pos = [0, -5];
obs_size = 1;
obs_width = 0.1;
obs_interval = [-45, 45];
obs = build_obs(obs_pos, obs_width, obs_size);

source_pos = [0, -5];
source_support_width = 10;
source_support_size = 2;
source_function = @(x_axis) (0.5*(x_axis < 0) + (1*x_axis >= 0));
%source_function = @(x_axis) (1);
source = build_source(source_pos, source_support_width, source_support_size, source_function);

[mirr_x, mirr_h, mirr_g] = rendering(obs, obs_interval, source, mirror_BRDF, step_lin, margin, false, step_ang, true, mirror_BRDF, margin, false);
[x, h, g] = rendering(obs, obs_interval, source, blurred_mirror_BRDF, step_lin, sigma, false, step_ang, true, mirror_BRDF, margin, false);

f_est1 = derendering(g, h, 1);
% f_est2 = derendering(g, h, 2);
% f_est3 = derendering(g, h, 3);
% f_est4 = derendering(g, h, 4);

figure;
subplot(1, 2, 1)
stem(mirr_x, mirr_g(length(mirr_g)/4:3*length(mirr_g)/4));
subplot(1, 2, 2)
stem(x, f_est1(length(f_est1)/4:3*length(f_est1)/4));

[f_est1, h_est1] = blind_derendering(g, 1000, 1000, 3);
%[f_est2, h_est2] = blind_derendering(g, 1000, 1000, 2);
%[f_est3, h_est3] = blind_derendering(g, 1000, 1000, 3);
%[f_est4, h_est4] = blind_derendering(g, 1000, 1000, 4);

figure;
subplot(1, 2, 1)
%stem(x, f_est1(length(f_est1)/4:3*length(f_est1)/4));
stem(f_est1);
subplot(1, 2, 2)
stem(h_est1);

% figure;
% subplot(1, 4, 1)
% stem(mirr_x, f_est1);
% subplot(1, 4, 2)
% stem(mirr_x, f_est2(length(f_est2)/4:3*length(f_est2)/4));
% subplot(1, 4, 3)
% stem(mirr_x, f_est3);
% subplot(1, 4, 4)
% stem(mirr_x, f_est4(length(f_est4)/4:3*length(f_est4)/4));
% 
% figure;
% subplot(1, 4, 1)
% stem(mirr_x, h_est1);
% subplot(1, 4, 2)
% stem(mirr_x, h_est2);
% subplot(1, 4, 3)
% stem(mirr_x, h_est3);
% subplot(1, 4, 4)
% stem(mirr_x, h_est4(length(h_est4)/4:3*length(h_est4)/4));

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

function obs = build_obs(obs_pos, obs_width, obs_size)
    obs = zeros(obs_size, 2);
    if(obs_size == 1)
        x_axis = zeros(1, 1);
    else
        x_axis = linspace(-obs_width/2, obs_width/2, obs_size);
    end
    for i=1:obs_size
        obs(i, :) = [obs_pos(1)+x_axis(i), obs_pos(2)];
    end
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
function [x_axis, brdf_map_wall, l_out] = rendering(obs, obs_interval, source, brdf, step_lin, param, test_conv1, step_ang, test_conv2, mirror_BRDF, margin, output_conv)

    l_out = zeros(size(obs, 1), round(1/step_lin));
    brdf_map = zeros(round(1/step_ang), round(1/step_lin));
    x_axis = zeros(size(obs, 1), round(1/step_lin));
    wall_points_global = zeros(size(obs, 1), 2, 2);    
    obs_angle = linspace(obs_interval(1)+90, obs_interval(2)+90, 2);
    brdf_map_wall = zeros(1, round(1/step_lin));
    
    for i=1:size(obs, 1)
        
        wall_points = zeros(2, 2);
        for j=1:2
            obs_direction = [cos(obs_angle(j)*pi/180), sin(obs_angle(j)*pi/180)];
            wall_points(j, :) = obs(i, :) + obs_direction*ray_wall_intersection(obs(i, :), obs_direction);
        end         
    
%         yline(0); hold on
%         scatter(wall_points(:, 1), wall_points(:, 2)); hold on
%         scatter(obs(i, 1), obs(i, 2)); hold on
%         for j=1:2
%             line([obs(i, 1), wall_points(j, 1)], [obs(i, 2), wall_points(j, 2)]); hold on
%         end        
        
        x_axis(i, :) = linspace(wall_points(1, 1), wall_points(2, 1), round(1/step_lin));
        l_out_temp = zeros(1, round(1/step_lin));
        
        if test_conv1
            c_test = zeros(1, round(1/step_lin));
            c_test2 = zeros(1, round(1/step_lin));
        end
        
        for j=1:round(1/step_lin) 
            normal = [x_axis(j), -1] - [x_axis(j), 0];
            dir_out = obs(i, :) - [x_axis(j), 0];
            angle_out = atan2d(normal(1)*dir_out(2)-normal(2)*dir_out(1),normal(1)*dir_out(1)+normal(2)*dir_out(2));
            %angle_out = atand(x_axis(j)/5);
            
            if test_conv1
                c_angles_storage = zeros(1,size(source, 1));
                c_test2(j) = angle_out;
            end
            
            for k=1:size(source, 1)                
                dir_in = [source(k, 1), source(k, 2)] - [x_axis(j), 0];
                angle_in = atan2d(normal(1)*dir_in(2)-normal(2)*dir_in(1),normal(1)*dir_in(1)+normal(2)*dir_in(2));  
                %angle_in = atand((source(k, 1) + x_axis(j))/5);
                l_out_temp(j) = l_out_temp(j) + source(k, 3)*brdf(angle_in + angle_out, param); 

                if test_conv1
                    c_angles_storage(k) = angle_in;
                end
                
            end
            
            angles = linspace(-90, 90, round(1/step_ang));
            for k = 1:length(angles)
                brdf_map(k, j) = brdf(-(angles(k) + angle_out), param);
            end   
            
            if test_conv1
                c_angles = linspace(-90, 90, round(1/step_ang));   
                c_lumin = zeros(1, length(c_angles));
                c_brdf = zeros(1, length(c_angles));
                for k = 1:length(c_angles)
                    c_brdf(k) = brdf_map(k, j);
                end
                for k=1:length(c_angles_storage)
                    [~,idx] = min(abs(c_angles-c_angles_storage(k)));
                    c_lumin(idx) = source(k, 3);
                end
                temp = conv(c_lumin, c_brdf);
                c_angles2 = linspace(-90, 90, length(temp));
                [~,idx] = min(abs(c_angles2+angle_out));
                c_test(end-j+1) = temp(idx);
            end           
            
        end 
        
        l_out(i, :) = l_out_temp;        
        wall_points_global(i, :, :) = wall_points;       
    end
    
    if test_conv2        
        l_out2 = zeros(size(obs, 1), round(1/step_lin));
        l_in_test = zeros(size(obs, 1), round(1/step_lin));
        
        for i=1:size(obs, 1)
            
            % CREATE INCOMING LIGHT MAP (SIMULATE MIRROR CASE)
            wall_points = zeros(2, 2);
            for j=1:2
                obs_direction = [cos(obs_angle(j)*pi/180), sin(obs_angle(j)*pi/180)];
                wall_points(j, :) = obs(i, :) + obs_direction*ray_wall_intersection(obs(i, :), obs_direction);
            end    
            x_axis(i, :) = linspace(wall_points(1, 1), wall_points(2, 1), round(1/step_lin));
            l_in_test_temp = zeros(1, round(1/step_lin));
            for j=1:round(1/step_lin)
                normal = [x_axis(j), -1] - [x_axis(j), 0];
                dir_out = obs(i, :) - [x_axis(j), 0];
                angle_out = atan2d(normal(1)*dir_out(2)-normal(2)*dir_out(1),normal(1)*dir_out(1)+normal(2)*dir_out(2));
                for k=1:size(source, 1)                
                    dir_in = [source(k, 1), source(k, 2)] - [x_axis(j), 0];
                    angle_in = atan2d(normal(1)*dir_in(2)-normal(2)*dir_in(1),normal(1)*dir_in(1)+normal(2)*dir_in(2));  
                    l_in_test_temp(j) = l_in_test_temp(j) + source(k, 3)*mirror_BRDF(angle_in+angle_out, margin);
                end
            end
            l_in_test(i, :) = l_in_test_temp; 
            
            % CREATE BRDF MAP
            angle_diff_vec = linspace(obs_interval(1), obs_interval(2), round(1/step_ang));
            brdf_map = brdf(angle_diff_vec, param); 
            brdf_map_wall = zeros(1, round(1/step_lin));
            for j=1:round(1/step_ang)   
                ang = angle_diff_vec(j);
                wp = obs(2)*tand(ang);
                [~,idx] = min(abs(x_axis-wp));
                brdf_map_wall(idx) = brdf_map(j);
            end            
            
            % FINAL CONVOLUTION
            l_out2 = conv(l_in_test(i, :), brdf_map_wall, "full"); 
            %l_out2(l_out2 > 1) = 0;
            
        end
    end
    
    % PLOTS
    
    if test_conv1
        figure;
        
        subplot(1,3,1)
        stem(x_axis, l_out);
        title("Observed signal (rendering)")
        
        subplot(1,3,2)
        imagesc(brdf_map);
        xlabel("Output Angle");
        ylabel("Incident Angle (from -90 to 90)");
        
        subplot(1,3,3)
        stem(abs(obs(2))*tand(c_test2) + obs(1), c_test);
        title("Observed signal (convolution)")
        fprintf("Convolution error : %d \n", immse(l_out, flip(c_test)));        
    end
       
    if test_conv2
        figure;
        
        subplot(1,3,1)
        stem(x_axis, l_in_test);
        title("Observed signal (mirror)")
        
        subplot(1,3,2)
        stem(x_axis, brdf_map_wall);
        xlabel("Angle Difference");
        ylabel("BRDF");
        
        subplot(1,3,3)
        stem(x_axis, l_out2(length(l_out2)/4:3*length(l_out2)/4));
        title("Observed signal (convolution)")
        
        figure;
        
        subplot(1,2,1)
        stem(x_axis, l_out);
        title("Observed signal (rendering)")
        
        subplot(1,2,2)
        stem(x_axis, l_out2(length(l_out2)/4:3*length(l_out2)/4));
        title("Observed signal (convolution)")
    end  
    
    if output_conv
        l_out = l_out2;
    else
        l_out = [zeros(1, length(l_out)/2) l_out zeros(1, length(l_out)/2 - 1)];
    end
end

function f_est = derendering(g, h, method)
     if method == 1
         f_est = deconvlucy(g, h, 10, 0.00001);
     elseif method == 2
         f_est = deconvwnr(g, h, 0.001);
     elseif method == 3
         f_est = deconvreg(g, h);
     else
         f_est = TikDeconv(g, h ,10);
     end
end

function [f_est, h_est] = blind_derendering(g, f_size, h_size, method)
    if method == 1
        [f_est, h_est] = FastMotionDeblurring(g, f_size, h_size, 20, 1, 1, 5, 0.5, 1, 5, 0.1, 0.1);
        [f_est, h_est] = LucyRichardsonBlind(g, f_est, h_est, 100, 10);
    elseif method == 2
        [f_est, h_est] = NormSparsBlindDeconv(g, f_size, h_size, 100, 20, 20, 0.5, 2, 2, 0.001, 1, 20, 10);
        [f_est, h_est] = LucyRichardsonBlind(g, f_est, h_est, 100, 10);
    elseif method == 3        
        [f_est, h_est] = IterativeBlindDeconv(g, h_size, f_size, 0.001, 5);
        [f_est, h_est] = LucyRichardsonBlind(g, f_est, h_est, 100, 10);
    elseif method == 4
        [f_est, h_est] = LucyRichardsonBlind(g, f_size, h_size, 100, 10);
    else
        init = ones(1, h_size);
        %init(round(length(init)/2)) = 1;
        [f_est, h_est] = deconvblind(g, init, 50);
    end
end





