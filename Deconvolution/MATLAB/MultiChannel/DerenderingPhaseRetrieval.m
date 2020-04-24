function [s_est, g_est] = DerenderingPhaseRetrieval(d_interf, s_est, g_est, T, n, tau, obs, empty_source, num_ang, tol)  
    R1 = inf;
    R2 = inf;
    deltaR = inf;
    while deltaR > tol
        % 1. Update s
        % 1.1 Optimize W with s
        sa = optimvar('sa', T);
        Rfun = @(sa) computeR(sa, g_est, d_interf, n, obs, num_ang, empty_source);
        Rexp = fcn2optimexpr(Rfun,sa);
        Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp);
        Rprob.Constraints.cons1 = sa(1:floor(T/2)) == flip(sa(ceil(T/2) + 1:end));
        %Rprob.Constraints.cons2 = sa(:) >= 0;
        %Rprob.Constraints.cons3 = sa(:) <= 1;
        R0.sa = s_est;            
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'Display','iter', 'MaxFunctionEvaluations', 1e5));
        % 1.2 Assignements
        s_est = abs(Rsol.sa);
        R1p = R1;
        R1 = Rfval;
        % 2. Update g
        % 2.1 Optimize W with g
        gij = optimvar('gij', tau, n);
        Rfun = @(gij) computeR(s_est, gij, d_interf, n, obs, num_ang, empty_source);
        Rexp = fcn2optimexpr(Rfun,gij);
        Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp); 
        %Rprob.Constraints.cons1 = gij(:) >= 0;
        %Rprob.Constraints.cons2 = gij(:) <= 1;
        R0.gij = g_est;            
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fminunc,'Display','iter', 'MaxFunctionEvaluations', 1e5));
        % 2.2 Assignements
        g_est = Rsol.gij;
        R2p = R2;
        R2 = Rfval;
        deltaR = max([R1p - R1 R2p - R2]);            
        disp(deltaR);
    end
end

function out = computeR(s_est, g_est, d_interf, n, obs, num_ang, empty_source)
    % 1.1.1 Compute V
    R = 0;
    for k=1:n
        for l=1:n
            g_est_interf = xcorr(g_est);
            s_est_interf = xcorr(s_est);
            temp_diff = d_interf(:, get_col_num(k,l,n)) - FastRenderingInterf(obs, g_est_interf(:, get_col_num(k,l,n)), s_est_interf, num_ang, empty_source);
            temp_diff = temp_diff .* temp_diff;
            R = R + sum(temp_diff);
        end
    end 
    out = R;
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end
