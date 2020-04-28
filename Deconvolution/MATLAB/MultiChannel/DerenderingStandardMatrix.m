function [s_est, g_est] = DerenderingStandardMatrix(d, s_est, g_est, T, n, tau, obs, empty_source, num_ang, tol)  
    R1 = inf;
    R2 = inf;
    deltaR = inf;
    alpha = 1;
    while deltaR > tol
        % 1. Update s
        % 1.1 Optimize W with s
        sa = optimvar('sa', T, tau);
        Rfun = @(sa) computeR(sa, g_est, d, obs, num_ang, empty_source, alpha);
        Rexp = fcn2optimexpr(Rfun,sa);
        Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp);
        %Rprob.Constraints.cons1 = sa(1:floor(T/2), :) == flip(sa(ceil(T/2)+1:end, :));
        Rprob.Constraints.cons2 = sa(:) >= 0;
        Rprob.Constraints.cons4 = sum(sa, 1)./sum(sa(:, 1), 1) == ones(1, tau);
        %Rprob.Constraints.cons3 = sa(:) <= 1;
        R0.sa = s_est;            
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 1e5));
        % 1.2 Assignements
        s_est = Rsol.sa;
        R1p = R1;
        R1 = Rfval;
        % 2. Update g
        % 2.1 Optimize W with g
        gij = optimvar('gij', tau, n);
        Rfun = @(gij) computeR(s_est, gij, d, obs, num_ang, empty_source, alpha);
        Rexp = fcn2optimexpr(Rfun,gij);
        Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp); 
        Rprob.Constraints.cons1 = gij(:) >= 0;
        %Rprob.Constraints.cons2 = gij(:) <= 1;
        R0.gij = g_est;            
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 1e5));
        % 2.2 Assignements
        g_est = Rsol.gij;
        R2p = R2;
        R2 = Rfval;
        deltaR = max([R1p - R1 R2p - R2]);            
        disp(deltaR);
    end
end

function out = computeR(s_est, g_est, d, obs, num_ang, empty_source, alpha)
    % 1.1.1 Compute V
    R = 0;
    for i=1:obs(5)
        temp = d(:, i) - s_est*g_est(:, i);
        temp_diff = temp .* temp;
        R = R + sum(temp_diff);
    end
    R = R + alpha*norm(s_est, 1);
    out = R;
end