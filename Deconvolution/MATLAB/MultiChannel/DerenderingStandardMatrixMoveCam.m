function [s_est, g_est] = DerenderingStandardMatrixMoveCam(d, s_est, g_est_summed, g_est, T, n, tau, obs, empty_source, num_ang, tol)
    R1 = inf;
    R2 = inf;
    deltaR1 = inf;
    deltaR2 = inf;
    alpha = 0;
    while deltaR1 + deltaR2 > tol
        % 1. Update s
        % 1.1 Optimize W with s
        for i=1:n
            sa = optimvar('sa', tau);
            Rfun = @(sa) computeR(sa, g_est(:, :, i), d(:, i), obs, num_ang, empty_source, alpha);
            Rexp = fcn2optimexpr(Rfun,sa);
            Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp);
            %Rprob.Constraints.cons1 = sa(1:floor(T/2), :) == flip(sa(ceil(T/2)+1:end, :));
            Rprob.Constraints.cons2 = sa(:) >= 0;
            %Rprob.Constraints.cons4 = sum(sa, 1)./sum(sa(:, 1), 1) == ones(1, tau);
            %Rprob.Constraints.cons3 = sa(:) <= 1;
            R0.sa = s_est;            
            [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 1e5, 'Display', 'off'));
            % 1.2 Assignements
            s_est = Rsol.sa;
            R1p = R1;
            R1 = Rfval;
            deltaR1 = abs(R1p - R1);            
        end
        % 2. Update g
        % 2.1 Optimize W with g
        for i=1:n
            gij = optimvar('gij', T, tau);
            Rfun = @(gij) computeR(s_est, gij, d(:, i), obs, num_ang, empty_source, alpha);
            Rexp = fcn2optimexpr(Rfun,gij);
            Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp); 
            Rprob.Constraints.cons1 = gij(:) >= 0;
            Rprob.Constraints.cons2 = sum(gij, 2) == g_est_summed(:, i);
            R0.gij = g_est(:, :, i);            
            [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 1e5, 'Display', 'off'));
            % 2.2 Assignements
            g_est(:, :, i) = Rsol.gij;
            R2p = R2;
            R2 = Rfval;  
            deltaR2 = abs(R2p - R2);            
        end
        disp(deltaR1 + deltaR2);
    end
end

function out = computeR(s_est, g_est, d, obs, num_ang, empty_source, alpha)
    % 1.1.1 Compute V
    R = 0;
    temp = d - g_est*s_est;
    temp_diff = temp .* temp;
    R = R + sum(temp_diff);
    R = R + alpha*norm(g_est, 1);
    out = R;
end