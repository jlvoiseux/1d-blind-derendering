function l_out = FastRenderingInterf(obs, source_interf, brdf_interf, num_ang, empty_source)
    T = length(brdf_interf);
    l_out = zeros(obs(5), T);
    angles_brdf_map = linspace(-180, 180, num_ang);
    x_axis = zeros(obs(5), num_lin);
    
    obs_angle = linspace(obs(3)+90, obs(4)+90, 2);
    wall_points = zeros(2, 2);    
    for j=1:2
        obs_direction = [cos(obs_angle(j)*pi/180), sin(obs_angle(j)*pi/180)];
        wall_points(j, :) = obs([1 2]) + obs_direction*ray_wall_intersection(obs([1 2]), obs_direction);
    end    
    for i=1:obs(5)
        x_axis(i, :) = linspace(wall_points(1, 1), wall_points(2, 1), num_lin);
        l_out_temp = zeros(1, T);
        for j=1:T
            normal = [x_axis(i, j), -1] - [x_axis(i, j), 0];
            dir_out = obs([1 2]) - [x_axis(i, j), 0];
            angle_out = atan2d(normal(1)*dir_out(2)-normal(2)*dir_out(1),normal(1)*dir_out(1)+normal(2)*dir_out(2));
            lumin_map = zeros(1, num_ang);
            for k=1:size(empty_source, 1)
                dir_in = [empty_source(k, 1, i), empty_source(k, 2, i)] - [x_axis(i, j), 0];
                angle_in = atan2d(normal(1)*dir_in(2)-normal(2)*dir_in(1),normal(1)*dir_in(1)+normal(2)*dir_in(2));
                [~,idx] = min(abs(angles_brdf_map-angle_in));
                lumin_map(idx) = lumin_map(idx) + empty_source(k, 3, i);
            end
            temp1 = xcorr(lumin_map);
            temp1(temp1 > 0.5) = source_interf(i);
            temp2 = xcorr(brdf_map);
            temp = conv(temp1, temp2);
            angles_temp = linspace(angles_brdf_map(1), angles_brdf_map(end), length(temp));
            [~,idx] = min(abs(angles_temp + angle_out));
            l_out_temp(end-j+1) = temp(idx);
        end
    end
    l_out(i, :) = l_out_temp;
end

