function [wall_points, out] = CreateRenderingMatrixFromBRDFFinal(obs, source, brdf, num_lin, param)
    
    k = size(source, 1); % Source Support Size
    T = num_lin;
    n = obs(5); % Number of moves    
    out = zeros(T, k, n);    

    wall_points = ComputeWallPoints(obs, T, n);    
    
    for i=1:n
        
        normal = [wall_points(i, :); -1*ones(1, T)] - [wall_points(i, :); zeros(1, T)];
        dir_out = [obs(1)*ones(1, T); obs(2)*ones(1, T)] - [wall_points(i, :); zeros(1, T)];
        angle_out = atan2d(normal(1, :).*dir_out(2, :)-normal(2, :).*dir_out(1, :), normal(1, :).*dir_out(1, :)+normal(2, :).*dir_out(2, :)); 
        
        for k=1:size(source, 1)                
            dir_in = [source(k, 1, 1).*ones(1, T) ; source(k, 2, 1).*ones(1, T)] - [wall_points(i, :) ; zeros(1, T)];
            angle_in = atan2d(normal(1, :).*dir_in(2, :)-normal(2, :).*dir_in(1, :), normal(1, :).*dir_in(1, :)+normal(2, :).*dir_in(2, :));  
            out(:, k, i) = brdf(angle_in, angle_out, param);                
        end
        
    end
    
end