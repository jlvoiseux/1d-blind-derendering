function out = ComputeWallPointsIds(wall_points)
    out = zeros(size(wall_points));
    nmove = size(wall_points, 1);
    T = size(wall_points, 2);
    init = round(nmove/2);
    for j=init:nmove
        wp = wall_points(j, :);
        if j ~= init
            wp_prev_large = repmat(prev_wp', 1, T);
            [~, out(j, :)] = min(abs(wp_prev_large - wp));
        end
        prev_wp = wp;
    end
    prev_wp = wall_points(init, :);
    for j=init-1:-1:1
        wp = wall_points(j, :);
        wp_prev_large = repmat(prev_wp', 1, T);
        [~, out(j, :)] = min(abs(wp_prev_large - wp));
        prev_wp = wp;
    end
end

