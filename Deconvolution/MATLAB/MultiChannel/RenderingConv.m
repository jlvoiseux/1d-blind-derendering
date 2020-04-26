%% Rendering through a single convolution of two quantities present on the wall

function l_out = RenderingConv(obs, source, brdf, param, num_lin, num_ang, mirror_BRDF, margin)
    x_axis = zeros(obs(5), num_lin);
    obs_angle = linspace(obs(3)+90, obs(4)+90, 2);
    l_out = zeros(obs(5), num_lin);
    l_in_test = zeros(obs(5), num_lin);
    angle_diff_vec = linspace(obs(3), obs(4), num_ang);
    brdf_map = brdf(angle_diff_vec, param); 
        
    for i=1:obs(5)
        % CREATE INCOMING LIGHT MAP (SIMULATE MIRROR CASE)
        wall_points = zeros(2, 2);
        for j=1:2
            obs_direction = [cos(obs_angle(j)*pi/180), sin(obs_angle(j)*pi/180)];
            wall_points(j, :) = obs([1 2]) + obs_direction*ray_wall_intersection(obs([1 2]), obs_direction);
        end    
        x_axis(i, :) = linspace(wall_points(1, 1), wall_points(2, 1), num_lin);
        l_in_test_temp = zeros(1, num_lin);
        for j=1:num_lin
            normal = [x_axis(i,j), -1] - [x_axis(i,j), 0];
            dir_out = obs([1 2]) - [x_axis(i,j), 0];
            angle_out = atan2d(normal(1)*dir_out(2)-normal(2)*dir_out(1),normal(1)*dir_out(1)+normal(2)*dir_out(2));
            for k=1:size(source, 1)                
                dir_in = [source(k, 1), source(k, 2)] - [x_axis(i,j), 0];
                angle_in = atan2d(normal(1)*dir_in(2)-normal(2)*dir_in(1),normal(1)*dir_in(1)+normal(2)*dir_in(2));  
                l_in_test_temp(j) = l_in_test_temp(j) + source(k, 3, i)*mirror_BRDF(angle_in+angle_out, margin);
            end
        end
        l_in_test(i, :) = l_in_test_temp;         
        
        % CREATE BRDF MAP        
        brdf_map_wall = zeros(1, num_lin);
        for j=1:num_ang  
            ang = angle_diff_vec(j);
            wp = obs(2)*tand(ang);
            [~,idx] = min(abs(x_axis(i, :)-wp));
            brdf_map_wall(idx) = brdf_map(j);
        end            
        % FINAL CONVOLUTION
        l_out(i, :) = conv(l_in_test(i, :), brdf_map_wall, "same");
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