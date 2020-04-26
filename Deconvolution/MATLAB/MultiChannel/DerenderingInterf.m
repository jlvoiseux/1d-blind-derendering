function [s_interf_est, g_interf_est] = DerenderingInterf(d_interf, s_interf_est, g_interf_est, T, n, tau, obs, empty_source, num_ang, tol)  
    R1 = inf;
    R2 = inf;
    deltaR = inf;
    alpha = [1e5];
    for i = 1:length(alpha)
        deltaR = inf;
        while deltaR > tol
            % 1. Update s
            % 1.1 Optimize W with s
            sa = optimvar('sa', 2*T-1);
            Rfun = @(sa) computeR(sa, g_interf_est, d_interf, obs, num_ang, tau, alpha(i), empty_source);
            Rexp = fcn2optimexpr(Rfun,sa);
            Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp);
            %Rprob.Constraints.cons1 = sa(T) == 1;            
            Rprob.Constraints.cons2 = sa(1:T-1) == flip(sa(T+1:end));
            %Rprob.Constraints.cons3 = sa(:) <= 1;
            %Rprob.Constraints.cons4 = sa(:) >= 0;
            R0.sa = s_interf_est;            
            [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'Display','iter', 'MaxFunctionEvaluations', 1e5));
            % 1.2 Assignements
            s_interf_est = Rsol.sa;
            s_interf_est = abs(s_interf_est);            
            R1p = R1;
            R1 = Rfval;
            % 2. Update g
            % 2.1 Optimize W with g
            gij = optimvar('gij', 2*tau-1, n*n);
            Rfun = @(gij) computeR(s_interf_est, gij, d_interf, obs, num_ang, tau, alpha(i), empty_source);
            Rexp = fcn2optimexpr(Rfun,gij);
            Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp); 
            %Rprob.Constraints.cons1 = gij(:) <= 1;  
            %Rprob.Constraints.cons2 = gij(:) >= 0;  
            R0.gij = g_interf_est;            
            [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fminunc,'Display','iter', 'MaxFunctionEvaluations', 1e5));
            % 2.2 Assignements
            g_interf_est = Rsol.gij;
            R2p = R2;
            R2 = Rfval;
            deltaR = R2p - R2;            
            disp(deltaR);
        end  
    end
    s_interf_est = s_interf_est./max(s_interf_est);
end

function out = computeR(s_interf_est, g_interf_est, d_interf, obs, num_ang, tau, alpha, empty_source)
    % 1.1.1 Compute V    
    temp1 = FastRenderingInterf(obs, g_interf_est, s_interf_est, num_ang, empty_source);
    temp_diff = d_interf - temp1;
    temp_diff = temp_diff .* temp_diff;
    R = sum(sum(temp_diff));
        
%     g_same = g_interf_est(:, get_col_num(1:n, 1:n, n));
%     t = (-tau+1:tau-1)';
%     Rterm = sum(sum(t.*t.*g_same.*g_same));
%     R = R + alpha*Rterm;    
    out = R;
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end