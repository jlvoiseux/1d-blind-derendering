function out = CreateRenderingMatrixFromBRDF(obs, source, brdf, num_lin, param)
    k=size(source, 1);
    T = num_lin;
    out = zeros(k, T);
    obs_angle = linspace(obs(3)+90, obs(4)+90, 2);
    wall_points = zeros(2, 2);    
    for j=1:2
        obs_direction = [cos(obs_angle(j)*pi/180), sin(obs_angle(j)*pi/180)];
        wall_points(j, :) = obs([1 2]) + obs_direction*ray_wall_intersection(obs([1 2]), obs_direction);
    end    
               
    x_axis = linspace(wall_points(1, 1), wall_points(2, 1), T); 
    normal = [x_axis; -1*ones(1, T)] - [x_axis; zeros(1, T)];
    dir_out = [obs(1)*ones(1, T); obs(2)*ones(1, T)] - [x_axis ; zeros(1, T)];
    angle_out = atan2d(normal(1, :).*dir_out(2, :)-normal(2, :).*dir_out(1, :), normal(1, :).*dir_out(1, :)+normal(2, :).*dir_out(2, :)); 
    for k=1:size(source, 1)                
        dir_in = [source(k, 1, 1).*ones(1, T) ; source(k, 2, 1).*ones(1, T)] - [x_axis ; zeros(1, T)];
        angle_in = atan2d(normal(1, :).*dir_in(2, :)-normal(2, :).*dir_in(1, :), normal(1, :).*dir_in(1, :)+normal(2, :).*dir_in(2, :));  
        out(k, :) = brdf(angle_in + angle_out, param);                
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
