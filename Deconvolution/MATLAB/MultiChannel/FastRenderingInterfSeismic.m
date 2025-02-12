%% Produces a normal render using normal values

function l_out = FastRenderingInterfSeismic(obs, source, brdf_interf, num_ang, empty_source)
    n = obs(5);
    T = length(brdf_interf);
    angles_brdf_map = linspace(-180, 180, num_ang).*ones(T, num_ang);
    x_axis = zeros(1, T);
    obs_angle = linspace(obs(3)+90, obs(4)+90, 2);
    wall_points = zeros(2, 2);    
    l_out = zeros(T, n*n);
    for j=1:2
        obs_direction = [cos(obs_angle(j)*pi/180), sin(obs_angle(j)*pi/180)];
        wall_points(j, :) = obs([1 2]) + obs_direction*ray_wall_intersection(obs([1 2]), obs_direction);
    end 
    
    x_axis(:) = linspace(wall_points(1, 1), wall_points(2, 1), T);
    l_out_temp = zeros(1, T);
    normal = [x_axis(1:T); -1*ones(1, T)] - [x_axis(1:T); zeros(1, T)];
    dir_out = [obs(1)*ones(1, T); obs(2)*ones(1, T)] - [x_axis(1:T); zeros(1, T)];
    angle_out = atan2d(normal(1, :).*dir_out(2, :)-normal(2, :).*dir_out(1, :), normal(1, :).*dir_out(1, :)+normal(2, :).*dir_out(2, :)); 
    lumin_map_obs = zeros(T, num_ang, n);
    for i=1:n
        lumin_map = zeros(T, num_ang);
        for k=1:size(empty_source, 1)
            dir_in = [empty_source(k, 1, 1)*ones(1, T); empty_source(k, 2, 1)*ones(1, T)] - [x_axis(1:T); zeros(1, T)];
            angle_in = atan2d(normal(1, :).*dir_in(2, :)-normal(2, :).*dir_in(1, :) ,normal(1, :).*dir_in(1, :)+normal(2, :).*dir_in(2, :));
            [~,idx] = min(abs(angles_brdf_map-angle_in'), [], 2);
            indexes = sub2ind(size(lumin_map),(1:T)',idx);
            lumin_map(indexes) = lumin_map(indexes) + source(k, i);
        end
        lumin_map_obs(:, :, i) = lumin_map;
    end
    temp = zeros(T, 1498, n*n);
    for j=1:T
        for i=1:n
            for k=1:n
                temp(j, :, get_col_num(i,k, n)) = conv(lumin_map_obs(j, :, k) + lumin_map_obs(j, :, i), brdf_interf); % Make full BRDF map to avoid 'same'
            end
        end
    end  
    angles_temp = linspace(angles_brdf_map(1, 1), angles_brdf_map(1, end), 1498).*ones(T, 1498);
    [~,idx] = min(abs(angles_temp + angle_out'), [], 2);
    indexes = sub2ind(size(angles_temp),(1:T)',idx);
    for i=1:n*n
        temp2 = temp(:, :, i);
        l_out_temp(end-(1:T)+1) = temp2(indexes);   
        l_out(:, i) = flip(l_out_temp)';
    end    
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
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