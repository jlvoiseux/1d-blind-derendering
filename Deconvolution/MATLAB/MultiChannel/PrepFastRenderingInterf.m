%% Produces an interferometric render using interferometric values

function [temp1, indexes] = PrepFastRenderingInterf(obs, num_ang, empty_source, tau, T)    
    angles_brdf_map = linspace(-180, 180, num_ang).*ones(2*T-1, num_ang);
    x_axis = zeros(1, 2*T-1);
    obs_angle = linspace(obs(3)+90, obs(4)+90, 2);
    wall_points = zeros(2, 2);      
    for j=1:2
        obs_direction = [cos(obs_angle(j)*pi/180), sin(obs_angle(j)*pi/180)];
        wall_points(j, :) = obs([1 2]) + obs_direction*ray_wall_intersection(obs([1 2]), obs_direction);
    end    
    
    x_axis(:) = linspace(wall_points(1, 1), wall_points(2, 1), 2*T-1);    
    normal = [x_axis(1:2*T-1); -1*ones(1, 2*T-1)] - [x_axis(1:2*T-1); zeros(1, 2*T-1)];
    dir_out = [obs(1)*ones(1, 2*T-1); obs(2)*ones(1, 2*T-1)] - [x_axis(1:2*T-1); zeros(1, 2*T-1)];
    angle_out = atan2d(normal(1, :).*dir_out(2, :)-normal(2, :).*dir_out(1, :), normal(1, :).*dir_out(1, :)+normal(2, :).*dir_out(2, :));
    lumin_map = zeros(2*T-1, num_ang);
    for k=1:size(empty_source, 1)
        dir_in = [empty_source(k, 1, 1)*ones(1, 2*T-1); empty_source(k, 2, 1)*ones(1, 2*T-1)] - [x_axis(1:2*T-1); zeros(1, 2*T-1)];
        angle_in = atan2d(normal(1, :).*dir_in(2, :)-normal(2, :).*dir_in(1, :) ,normal(1, :).*dir_in(1, :)+normal(2, :).*dir_in(2, :));
        [~,idx] = min(abs(angles_brdf_map-angle_in'), [], 2);
        indexes = sub2ind(size(lumin_map),(1:2*T-1)',idx);
        lumin_map(indexes) = lumin_map(indexes) + empty_source(k, 3, 1);
    end
    temp1 = zeros(2*T-1-T, 2*num_ang-1);    
    for j=round(T/2):2*T-1-round(T/2)
        temp1 = xcorr(lumin_map(j, :));
        temp1 = trim(temp1, 2*tau-1);        
    end
    angles_temp = linspace(angles_brdf_map(1, 1), angles_brdf_map(1, end), 2*(2*num_ang-1)-1).*ones(2*T-1, 2*(2*num_ang-1)-1);
    [~,indexes] = min(abs(angles_temp + angle_out'), [], 2);
    %indexes = sub2ind(size(angles_temp),(1:2*T-1)',idx);    
    
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

function out = trim(signal, num)
    [~, ind] = maxk(signal, num);
    out = zeros(size(signal));
    out(ind) = 1;
end

