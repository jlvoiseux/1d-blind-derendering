function g_interf_est = RecovExactXcorrMoveCam(initial_weights, d_interf, g_interf_est, nmove, T, tol)
    X1 = inf;
    X2 = inf;
    deltaX = inf;
    we = optimvar('we', nmove*nmove);
    xmat = optimvar('xmat', 2*T-1, nmove*nmove);
    Xfun = @(we, xmat) computeWeights(nmove, d_interf, xmat, we);
    Xexp = fcn2optimexpr(Xfun, we, xmat);
    Xprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Xexp);
    Xprob.Constraints.cons1 = we(:) >= min(initial_weights);
    X0.we = initial_weights;       
    X0.xmat = g_interf_est;
    [Xsol,Xfval,~,~] = solve(Xprob,X0,'Options', optimoptions(@fmincon,'Display','iter', 'MaxFunctionEvaluations', 1e5, 'MaxIterations', 1e6));
    initial_weights = Xsol.we;
    g_interf_est = Xsol.xmat;
    X1p = X1;
    X1 = Xfval;

%         we = optimvar('gi', 2*T-1, n*n);
%         Xfun = @(we) computeWeights(n, d_interf, we, initial_weights);
%         Xexp = fcn2optimexpr(Xfun,we);
%         Xprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Xexp);
%         X0.we = g_interf_est;            
%         [Xsol,Xfval,~,~] = solve(Xprob,X0,'Options', optimoptions(@fminunc,'Display','off', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e6));
%         g_interf_est = Xsol.we;
%         X2p = X2;
%         X2 = Xfval;
        

end

function out = computeWeights(n, d_interf, g_interf_est, weights)
    temp_diff = d_interf-weights'.*g_interf_est;
    temp_diff = temp_diff.*temp_diff;
    X = sum(sum(temp_diff));
    out = X;
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end