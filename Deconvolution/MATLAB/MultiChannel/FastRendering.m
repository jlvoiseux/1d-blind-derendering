function l_out = FastRendering(obs, source, brdf, num_ang, empty_source)
    T = length(brdf);
    angles_brdf_map = linspace(-180, 180, num_ang);
    x_axis = zeros(1, T);
    obs_angle = linspace(obs(3)+90, obs(4)+90, 2);
    wall_points = zeros(2, 2);    
    for j=1:2
        obs_direction = [cos(obs_angle(j)*pi/180), sin(obs_angle(j)*pi/180)];
        wall_points(j, :) = obs([1 2]) + obs_direction*ray_wall_intersection(obs([1 2]), obs_direction);
    end 
    
    x_axis(:) = linspace(wall_points(1, 1), wall_points(2, 1), T);
    l_out_temp = zeros(1, T);
    for j=1:T
        normal = [x_axis(j), -1] - [x_axis(j), 0];
        dir_out = obs([1 2]) - [x_axis(j), 0];
        angle_out = atan2d(normal(1)*dir_out(2)-normal(2)*dir_out(1),normal(1)*dir_out(1)+normal(2)*dir_out(2));
        lumin_map = zeros(1, num_ang);
        for k=1:size(empty_source, 1)
            dir_in = [empty_source(k, 1, 1), empty_source(k, 2, 1)] - [x_axis(j), 0];
            angle_in = atan2d(normal(1)*dir_in(2)-normal(2)*dir_in(1),normal(1)*dir_in(1)+normal(2)*dir_in(2));
            [~,idx] = min(abs(angles_brdf_map-angle_in));
            lumin_map(idx) = lumin_map(idx) + empty_source(k, 3, 1);
        end
        temp1 = lumin_map;
        temp1(temp1 > 0.5) = source;
        temp2 = brdf;
        temp = conv(temp1, temp2, 'same');
        angles_temp = linspace(angles_brdf_map(1), angles_brdf_map(end), length(temp));
        [~,idx] = min(abs(angles_temp + angle_out));
        l_out_temp(end-j+1) = temp(idx);
    end
    
    l_out = l_out_temp';
    
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