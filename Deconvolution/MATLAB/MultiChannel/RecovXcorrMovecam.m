function [weights, g_interf_est] = RecovXcorrMovecam(d, d_interf, nmove, nsource, T)
    g_interf_est_temp = zeros(2*T-1, nmove*nmove, nsource);
    weights = zeros(nsource, nmove*nmove);
    for h=1:nsource
        for i=1:nmove
            for j=1:nmove
                single_weight =  sqrt(mean(d(:, i, h))*mean(d(:, j, h)));
                weights(h, get_col_num(i,j, nmove)) = single_weight;
                g_interf_est_temp(:, get_col_num(i, j, nmove), nsource) = d_interf(:, get_col_num(i, j, nmove), h)/single_weight;
            end
        end
    end
    g_interf_est = reshape(mean(g_interf_est_temp, 3), [2*T-1, nmove*nmove]);
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end