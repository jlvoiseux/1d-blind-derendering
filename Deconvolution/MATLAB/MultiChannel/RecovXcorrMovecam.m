function g_interf_est = RecovXcorrMovecam(d, d_interf, n, T, tau)
    g_interf_est = zeros(2*T-1, n*n);
    for i=1:n
        for j=1:n
            g_interf_est(:, get_col_num(i, j, n)) = d_interf(:, get_col_num(i, j, n))/(mean(d(:, i))*mean(d(:, j)));
        end
    end
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end