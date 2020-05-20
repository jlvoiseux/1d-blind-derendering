function out = ComputeWallPoints(obs, T, n)
    
    out = zeros(n, T);
    
    % Compute start angles
    gap_start = obs(3) - obs(7);
    gap_end = obs(4);
    start_angles = linspace(gap_start, gap_end, n) + obs(7);
    sensor_length = 2*obs(8)*tand(obs(7)/2);
    
    for i=1:n
        angles = zeros(1, T);
        for j=1:T
            if j==1
                angles(j) = start_angles(i);
            else
                angles(j) = angles(1) - obs(7)/2 + atand(((sensor_length/2) - ((j-1)*sensor_length)/(T-1))/obs(8));
            end
            out(i, j) = obs(1)+obs(2)*tand(angles(j));
        end
    end
    
end