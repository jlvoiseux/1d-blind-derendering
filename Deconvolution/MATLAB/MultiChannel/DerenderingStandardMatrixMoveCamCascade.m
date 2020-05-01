function [s_est_full, g_est_deflat] = DerenderingStandardMatrixMoveCamCascade(d_full, s_est_full, g_est, T, nmove, nsource, tau, tol)
    g_est_flat = reshape(g_est, [T, tau*nmove]);
    d_full_flat = reshape(d_full, [T, nmove*nsource]);
    ac_mat = buildAutocorrMat(d_full, T, nsource, nmove);
    R1 = inf;
    R2 = inf;
    deltaR = inf;
    alpha = 0.1;
    while deltaR > tol
        % 1. Update s
        % 1.1 Optimize W with s
        sa = optimvar('sa', tau, nsource);
        Rfun = @(sa) computeR(sa, g_est_flat, d_full_flat, alpha, tau, nmove, nsource, ac_mat);
        Rexp = fcn2optimexpr(Rfun,sa);
        Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp);
        Rprob.Constraints.cons1 = sa(:) >= 0;
        Rprob.Constraints.cons2 = sa(:) <= 1;
        R0.sa = s_est_full;            
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 1e5, 'Display', 'iter'));
        % 1.2 Assignements
        s_est_full = Rsol.sa;
        R1p = R1;
        R1 = Rfval;
        for i=1:nmove-1
            % 2. Update g
            % 2.1 Optimize W with g
            Rfsum = 0;
            gij1 = optimvar('gij1', T, tau);
            gij2 = optimvar('gij2', T, tau);
            Rfun = @(gij1, gij2) computeRCascade(s_est_full, gij1, gij2, d_full_flat, tau, nmove, nsource, ac_mat(), i);
            Rexp = fcn2optimexpr(Rfun,gij1, gij2);
            Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp); 
            Rprob.Constraints.cons1 = gij1(:) >= 0;
            Rprob.Constraints.cons2 = gij1(:) <= 1;
            Rprob.Constraints.cons1 = gij2(:) >= 0;
            Rprob.Constraints.cons2 = gij2(:) <= 1;
            R0.gij1 = g_est_flat(:, 1+(i-1)*tau:i*tau); 
            R0.gij2 = g_est_flat(:, 1+(i)*tau:(i+1)*tau); 
            [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 1e5, 'Display', 'iter', 'UseParallel', false));
            % 2.2 Assignements
            g_est_flat(:, 1+(i-1)*tau:i*tau) = Rsol.gij1;
            g_est_flat(:, 1+(i)*tau:(i+1)*tau) = Rsol.gij2;
            Rfsum = Rfsum + Rfval;
        end
        R2p = R2;
        R2 = Rfsum;         
        deltaR = abs(max([R1p - R1 R2p - R2]));            
        disp(deltaR);
    end
    g_est_deflat = reshape(g_est_flat, [T, tau, nmove]);
end

function out = computeR(s_est_full, g_est_flat, d_full_flat, alpha, tau, nmove, nsource, ac_mat)
    % 1.1.1 Compute V
    R = 0;
    for i=1:nsource
        s_est = s_est_full(:, i);
        d = d_full_flat(:, 1+(i-1)*nmove:i*nmove);
        for j=1:nmove
            g_est = g_est_flat(:, 1+(j-1)*tau:j*tau);
            temp = d(:, j) - g_est*s_est;
            temp_diff = temp .* temp;
            R = R + sum(temp_diff);            
        end
    end
    out = R;
end

function out = computeRCascade(s_est_full, g_est1, g_est2, d_full_flat, tau, nmove, nsource, ac_mat, current_move)
    % 1.1.1 Compute V
    R = 0;
    for i=1:nsource
        s_est = s_est_full(:, i);
        d = d_full_flat(:, 1+(i-1)*nmove:i*nmove);
        for j=1:2
            ind = current_move-1+j;
            if j==1
                temp = d(:, ind) - g_est1*s_est;
                temp_diff = temp .* temp;
                R = R + sum(temp_diff);
            elseif j == 2
                temp = d(:, ind) - g_est2*s_est;
                temp_diff = temp .* temp;
                R = R + sum(temp_diff);
                num = max(ac_mat(ind-1, i), 1);
                for k=1:tau
                    R = R + sum(abs(g_est2(1:end-num+1, k)-g_est1(num:end)));
                end
            end            
        end
    end
    out = R;
end

function out = buildAutocorrMat(d_full, T, nsource, nmove)
    out = zeros(nmove-1, nsource);
    for i=1:nsource
        for j=2:nmove
            [~, ind] = max(xcorr(d_full(:, j-1, i), d_full(:, j, i)));
            out(j-1, i) = ind-T;
        end
    end
end