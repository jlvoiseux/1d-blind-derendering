function [s_interf_est, g_interf_est] = DerenderingInterf(d_interf, s_interf_est, g_interf_est, T, n, tau, obs, empty_source, num_ang, tol)  
    R1 = inf;
    R2 = inf;
    deltaR = inf;
    while deltaR > tol
        % 1. Update s
        % 1.1 Optimize W with s
        sa = optimvar('sa', 2*T-1);
        Rfun = @(sa) computeR(sa, g_interf_est, d_interf, n, obs, num_ang, empty_source);
        Rexp = fcn2optimexpr(Rfun,sa);
        Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp);
        Rprob.Constraints.cons1 = sa(T) == 1;            
        Rprob.Constraints.cons2 = sa(1:T-1) == flip(sa(T+1:end));
        Rprob.Constraints.cons3 = sa(:) <= 1;
        Rprob.Constraints.cons4 = sa(:) >= 0;
        R0.sa = s_interf_est;            
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'Display','iter', 'MaxFunctionEvaluations', 1e5));
        % 1.2 Assignements
        s_interf_est = Rsol.sa;
        R1p = R1;
        R1 = Rfval;
        % 2. Update g
        % 2.1 Optimize W with g
%         gij = optimvar('gij', 2*tau-1, n*n);
%         Rfun = @(gij) computeR(s_interf_est, gij, d_interf, n, obs, num_ang, empty_source);
%         Rexp = fcn2optimexpr(Rfun,gij);
%         Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp); 
%         Rprob.Constraints.cons1 = gij(:) <= 1;  
%         Rprob.Constraints.cons2 = gij(:) >= 0.01;  
%         R0.gij = g_interf_est;            
%         [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'Display','iter', 'MaxFunctionEvaluations', 1e5));
%         % 2.2 Assignements
%         g_interf_est = Rsol.gij;
%         R2p = R2;
%         R2 = Rfval;
%         deltaR = R2p - R2;            
%         disp(deltaR);
    end
end

function out = computeR(s_interf_est, g_interf_est, d_interf, n, obs, num_ang, empty_source)
    % 1.1.1 Compute V
    R = 0;
    for k=1:n
        for l=k:n
            temp_diff = d_interf(:, get_col_num(k, l, n)) - FastRenderingInterf(obs, g_interf_est(:, get_col_num(k, l, n)), s_interf_est, num_ang, empty_source);
            temp_diff = temp_diff .* temp_diff;
            R = R + sum(temp_diff);
        end
    end 
    out = R;
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end