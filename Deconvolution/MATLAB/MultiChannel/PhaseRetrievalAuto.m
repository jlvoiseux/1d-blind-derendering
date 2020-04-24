function s_est = PhaseRetrievalAuto(s_interf_est, T)
    s_est = rand(T, 1);
    sa = optimvar('sa', T);
    Xfun = @(sa) computeX(s_interf_est, sa);
    Xexp = fcn2optimexpr(Xfun,sa);
    Xprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Xexp);
    Xprob.Constraints.cons1 = sa(:) >= 0;
    Xprob.Constraints.cons2 = sa(1:floor(T/2)) == flip(sa(ceil(T/2) + 1:end));
    Xprob.Constraints.cons3 = sa(:) <= 1;
    X0.sa = s_est;            
    [Xsol,~,~,~] = solve(Xprob,X0,'Options', optimoptions(@fmincon,'Display','iter', 'MaxFunctionEvaluations', 1e6, 'MaxIterations', 1e6));
    s_est = Xsol.sa;
end

function out = computeX(s_interf_est, s_est)
   temp_diff = s_interf_est - xcorr(s_est);
   temp_diff = temp_diff .* temp_diff;
   X = sum(temp_diff); 
   out = X;
end
