function [g_cov, s_est_full] = EstimateSourceMoveCam(d_full, s_est_full, g_est, T, nmove, nsource, tau, tol)
    d_cov = zeros(nmove, nsource);
    g_cov = zeros(tau, tau*nmove);
    for i=1:nmove
        g_cov(:, 1+(i-1)*tau:i*tau) = rand(tau, tau);
        for j=1:nsource
            d_cov(i, j) = cov(d_full(:, i, j));
        end
    end
    [dg, ndg] = buildindexlist(g_cov, tau, nmove);
    R1 = inf;
    R2 = inf;
    deltaR = inf;
    alpha = [0];
    for i=1:length(alpha)
        while deltaR > tol
            % 1. Update s
            % 1.1 Optimize W with s
            sa = optimvar('sa', tau, nsource);
            Rfun = @(sa) computeR(sa, g_cov, d_cov, alpha(i), tau, nmove, nsource);
            Rexp = fcn2optimexpr(Rfun,sa);
            Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp);
            Rprob.Constraints.cons2 = sa(:) >= 0;
            Rprob.Constraints.cons3 = sa(:) <= 1;
            R0.sa = s_est_full;            
            [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 1e5, 'Display', 'off'));
            % 1.2 Assignements
            s_est_full = Rsol.sa;
            R1p = R1;
            R1 = Rfval;
            % 2. Update g
            % 2.1 Optimize W with g
            gij = optimvar('gij', tau, tau*nmove);
            Rfun = @(gij) computeR(s_est_full, gij, d_cov, alpha(i), tau, nmove, nsource);
            Rexp = fcn2optimexpr(Rfun,gij);
            Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp); 
            Rprob.Constraints.cons1 = gij(ind2sub(size(gij),ndg)) >= -0.01;
            Rprob.Constraints.cons2 = gij(ind2sub(size(gij),ndg)) <= 0;
            Rprob.Constraints.cons3 = gij(ind2sub(size(gij),dg)) >= 0;
            Rprob.Constraints.cons4 = gij(ind2sub(size(gij),dg)) <= 0.1;
            %Rprob.Constraints.cons3 = gij(:) <= 0.1;
            %Rprob.Constraints.cons4 = gij(:) >= -0.01;
            R0.gij = g_cov;            
            [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 1e5, 'Display', 'off'));

            % 2.2 Assignements
            g_cov = Rsol.gij;
            R2p = R2;
            R2 = Rfval;         
            deltaR = abs(max([R1p - R1 R2p - R2]));            
            disp(deltaR);
        end
    end
end

function out = computeR(s_est_full, g_cov, d_cov, alpha, tau, nmove, nsource, ndg)
    % 1.1.1 Compute V
    R = 0;
    ndg = ~eye(tau);
    for i=1:nsource
        s_est = s_est_full(:, i);
        for j=1:nmove
            gc = g_cov(:, 1+(j-1)*tau:j*tau);
            temp_diff = d_cov(j, i) - s_est'*gc*s_est;
            temp_diff = temp_diff .* temp_diff;
            R = R + temp_diff; 
            R = R + alpha*abs(sum(sum(gc.*ndg)));
        end
    end
    %R = R + alpha*norm(g_cov, 1);
    
    out = R;
end

function [dg, ndg] = buildindexlist(g_cov, tau, nmove)
    idxdiag = zeros(tau, tau*nmove);
    idxnondiag = zeros(tau, tau*nmove);
    for j=1:nmove
        idxdiag(:, 1+(j-1)*tau:j*tau) = eye(tau);
        idxnondiag(:, 1+(j-1)*tau:j*tau) = ~eye(tau);
    end
    dg = find(idxdiag == 1);
    ndg = find(idxnondiag == 1);
end

function out = CheckCovMatrixEquality(g_cov, nmove, tau)
    for j=1:nmove
        gc = g_cov(:, 1+(j-1)*tau:j*tau);
        dgc = gc(logical(eye(tau)));
        dgcl = length(dgc(dgc <= 0));        
        if dgcl >= 1
            out = false;
            return
        end        
    end
    out = true;
end


