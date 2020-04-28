function s_est = PhaseRetrievalAuto(s_interf_est, T, tau)
    s_est = rand(T, tau);
    alpha = 0.01;
    sa = optimvar('sa', T, tau);
    Xfun = @(sa) computeX(s_interf_est, sa, tau, alpha);
    Xexp = fcn2optimexpr(Xfun,sa);
    Xprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Xexp);
    Xprob.Constraints.cons1 = sa(:) >= 0;
    Xprob.Constraints.cons2 = sum(sa, 1)./sum(sa(:, 1), 1) == ones(1, tau);
    %Xprob.Constraints.cons3 = sa(:) <= 1;
    X0.sa = s_est;            
    [Xsol,~,~,~] = solve(Xprob,X0,'Options', optimoptions(@fmincon,'Display','iter', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e6));
    s_est = Xsol.sa;
end

function out = computeX(s_interf_est, s_est, n, alpha)
    X = 0;
    for i=1:n
        temp_diff = s_interf_est(:, i) - xcorr(s_est(:, i));
        temp_diff = temp_diff .* temp_diff;
        X = X + sum(temp_diff);
    end   
    X = X + alpha*norm(s_est, 1);
   out = X;
end
