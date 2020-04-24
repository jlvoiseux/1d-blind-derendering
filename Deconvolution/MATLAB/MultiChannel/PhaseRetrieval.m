function g_est = PhaseRetrieval(g_interf_est, tau, n)
    g_est = rand(tau, n);
    gi = optimvar('gi', tau, n);
    Xfun = @(gi) computeX(n, g_interf_est, gi);
    Xexp = fcn2optimexpr(Xfun,gi);
    Xprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Xexp);
    Xprob.Constraints.cons1 = gi(:) >= 0;
    X0.gi = g_est;            
    [Xsol,~,~,~] = solve(Xprob,X0,'Options', optimoptions(@fmincon,'Display','iter', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e6));
    g_est = Xsol.gi;
end

function out = computeX(n, g_interf_est, g_est)
    X = 0;
    for k=1:n
        for l=1:n
           temp_diff = g_interf_est(:, get_col_num(k, l, n)) - xcorr(g_est(:,k), g_est(:,l));
           temp_diff = temp_diff .* temp_diff;
           X = X + sum(temp_diff);
        end
    end
    out = X;
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end