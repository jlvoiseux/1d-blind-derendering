function [s_est, g_est] = DerenderingPhaseRetrieval(d_interf, s_est, g_est, T, n, tau, obs, empty_source, num_ang, tol)  
    R1 = inf;
    R2 = inf;    
    alpha = [1e5, 0];
    for i=1:length(alpha)
        deltaR = inf;
        while deltaR > tol 
            % 1. Update s
            % 1.1 Optimize W with s
            sa = optimvar('sa', T);
            Rfun = @(sa) computeR(sa, g_est, d_interf, obs, num_ang, empty_source, tau, alpha(i), n);
            Rexp = fcn2optimexpr(Rfun,sa);
            Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp);
            Rprob.Constraints.cons1 = sa(1:floor(T/2)) == flip(sa(ceil(T/2) + 1:end));
            %Rprob.Constraints.cons2 = sa(:) >= 0;
            %Rprob.Constraints.cons3 = sa(:) <= 1;
            R0.sa = s_est;            
            [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'Display','iter', 'MaxFunctionEvaluations', 1e5));
            % 1.2 Assignements
            s_est = Rsol.sa;
            R1p = R1;
            R1 = Rfval;
            % 2. Update g
            % 2.1 Optimize W with g
            gij = optimvar('gij', tau, n);
            Rfun = @(gij) computeR(s_est, gij, d_interf, obs, num_ang, empty_source, tau, alpha(i), n);
            Rexp = fcn2optimexpr(Rfun,gij);
            Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp); 
            Rprob.Constraints.cons1 = gij(:) >= 0.1;
            Rprob.Constraints.cons2 = gij(:) <= 1;
            R0.gij = g_est;            
            [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'Display','iter', 'MaxFunctionEvaluations', 1e5));
            % 2.2 Assignements
            g_est = Rsol.gij;
            R2p = R2;
            R2 = Rfval;
            deltaR = max([R1p - R1 R2p - R2]);            
            disp(deltaR);
        end
    end
end

function out = computeR(s_est, g_est, d_interf, obs, num_ang, empty_source, tau, alpha, n)
    % 1.1.1 Compute V
    temp1 = FastRendering(obs, g_est, s_est, num_ang, empty_source);
    temp_diff = d_interf - xcorr(temp1);
    temp_diff = temp_diff .* temp_diff;
    R = sum(sum(temp_diff));

    g_interf_est = xcorr(g_est);
    g_same = g_interf_est(:, get_col_num(1:n, 1:n, n));
    t = (-tau+1:tau-1)';
    Rterm = sum(sum(t.*t.*g_same.*g_same));
    R = R + alpha*Rterm; 
    out = R;
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end