function g_est = PhaseRetrievalSpecial(g_interf_est, tau, n, sym)
    g_est = rand(tau, n);
    gi = optimvar('gi', tau, n);
    Xfun = @(gi) computeXspecial(n, g_interf_est, gi);
    Xexp = fcn2optimexpr(Xfun,gi);
    Xprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Xexp);
    Xprob.Constraints.cons1 = gi(:) >= 0;
    if sym
        Xprob.Constraints.cons2 = gi(1:floor(tau/2), :) == flip(gi(ceil(tau/2)+1:end, :));
    end
    X0.gi = g_est;            
    [Xsol,~,~,~] = solve(Xprob,X0,'Options', optimoptions(@fmincon,'Display','iter', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e6));
    g_est = Xsol.gi;
end

function out = computeXspecial(n, g_interf_est, g_est)
    X = 0;
    temp1 = xcorr(g_est);
    temp2 = sum(temp1, 2);
    temp3 = sum(g_interf_est - temp2);
    out = temp3;
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end