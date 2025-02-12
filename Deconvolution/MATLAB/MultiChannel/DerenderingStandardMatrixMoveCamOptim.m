function [s_est_full, g_est_deflat, g_est_flat] = DerenderingStandardMatrixMoveCamOptim(d_full, s_est_full, g_est, T, nmove, nsource, tau, tol)
    g_est_flat = reshape(g_est, [T, tau*nmove]);
    d_full_flat = reshape(d_full, [T, nmove*nsource]);
    ac_mat = buildAutocorrMat(d_full, T, nsource, nmove);
    ac_mat = round(mean(ac_mat, 2));
    R1 = inf;
    R2 = inf;
    deltaR = inf;
    alpha = 1;
    g_est_opt = buildOptFromFlat(g_est_flat, T, nmove, tau, ac_mat);
    while deltaR > tol
        % 1. Update g
        % 1.1 Optimize W with g         
        gij = optimvar('gij', size(g_est_opt, 1), tau);
        Rfun = @(gij) computeROpt(s_est_full, gij, d_full_flat, alpha, tau, nmove, nsource, ac_mat, T);
        Rexp = fcn2optimexpr(Rfun,gij);
        Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp); 
        Rprob.Constraints.cons1 = gij(:) >= 0;
        Rprob.Constraints.cons2 = gij(:) <= 1;
        R0.gij = g_est_opt;            
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 2e5, 'Display', 'iter', 'UseParallel', false, 'Algorithm', 'sqp'));
        % 2.2 Assignements
        g_est_opt = Rsol.gij;
        R2p = R2;
        R2 = Rfval;   
        % 2. Update s
        % 2.1 Optimize W with s
        sa = optimvar('sa', tau, nsource);
        Rfun = @(sa) computeROpt(sa, g_est_opt, d_full_flat, alpha, tau, nmove, nsource, ac_mat, T);
        Rexp = fcn2optimexpr(Rfun,sa);
        Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp);
        Rprob.Constraints.cons1 = sa(:) >= 0;
        Rprob.Constraints.cons2 = sa(:) <= 1;
        R0.sa = s_est_full;            
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 1e5, 'Display', 'iter', 'Algorithm', 'sqp'));
        % 1.2 Assignements
        s_est_full = Rsol.sa;
        R1p = R1;
        R1 = Rfval;              
        deltaR = max([R1p - R1 R2p - R2]);            
        disp(deltaR);
    end
    deltaR = inf;
    R1 = inf;
    R2 = inf;
    g_est_flat = buildFlatFromOpt(g_est_opt, T, nmove, tau, ac_mat);
    while deltaR > tol        
        % 1. Update g
        % 1.1 Optimize W with g
        gij = optimvar('gij', T, tau*nmove);
        Rfun = @(gij) computeR(s_est_full, gij, d_full_flat, alpha, tau, nmove, nsource);
        Rexp = fcn2optimexpr(Rfun,gij);
        Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp); 
        Rprob.Constraints.cons1 = gij(:) >= 0;
        Rprob.Constraints.cons2 = gij(:) <= 1;
        R0.gij = g_est_flat;            
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 2e5, 'Display', 'iter', 'UseParallel', true, 'Algorithm', 'sqp'));
        % 2.2 Assignements
        g_est_flat = Rsol.gij;
        R2p = R2;
        R2 = Rfval;       
        % 2. Update s
        % 2.1 Optimize W with s
        sa = optimvar('sa', tau, nsource);
        Rfun = @(sa) computeR(sa, g_est_flat, d_full_flat, alpha, tau, nmove, nsource);
        Rexp = fcn2optimexpr(Rfun,sa);
        Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp);
        Rprob.Constraints.cons1 = sa(:) >= 0;
        Rprob.Constraints.cons2 = sa(:) <= 1;
        R0.sa = s_est_full;            
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 1e5, 'Display', 'iter', 'Algorithm', 'sqp'));
        % 1.2 Assignements
        s_est_full = Rsol.sa;
        R1p = R1;
        R1 = Rfval;
        deltaR = max([R1p - R1 R2p - R2]);            
        disp(deltaR);
    end
    
    g_est_deflat = reshape(g_est_flat, [T, tau, nmove]);
end

function out = computeROpt(s_est_full, g_est_opt, d_full_flat, alpha, tau, nmove, nsource, ac_mat, T)
    % 1.1.1 Compute V
    R = 0;
    g_est_flat = buildFlatFromOpt(g_est_opt, T, nmove, tau, ac_mat);
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
    R = R + alpha*norm(g_est_opt, 1);
    out = R;
end

function out = computeR(s_est_full, g_est_flat, d_full_flat, alpha, tau, nmove, nsource)
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
    R = R + alpha*norm(g_est_flat, 1);
    out = R;
end

function out = buildAutocorrMat(d_full, T, nsource, nmove)
    out = zeros(nmove-1, nsource);
    for i=1:nsource
        for j=2:nmove
            [~, ind] = max(xcorr(d_full(:, j-1, i), d_full(:, j, i)));
            out(j-1, i) = max(0, ind-T);
        end
    end
end

function out = buildOptFromFlat(g_est_flat, T, nmove, tau, ac_mat)    
    length_util = sum(ac_mat > 0);
    out = zeros(T+sum(ac_mat)-length_util, tau);
    for j=1:nmove
        g_est = g_est_flat(:, 1+(j-1)*tau:j*tau);
        if j==1
            out(1:T, :) = g_est;
            curr_ind = T+1;
        else
            num = ac_mat(j-1);
            out(curr_ind:curr_ind+num, :) = g_est(end-num:end, :);
            curr_ind = curr_ind+num+1;
        end
    end
end

function out = buildFlatFromOpt(g_est_opt, T, nmove, tau, ac_mat)
    out = zeros(T, tau*nmove);
    for j=1:nmove        
        if j==1
            out(:, 1+(j-1)*tau:j*tau) = g_est_opt(1:T, :);
            g_prev = g_est_opt(1:T, :);
            curr_ind = T+1;
        else
            num = ac_mat(j-1);
            g_est = zeros(T, tau);
            g_est(1:end-num, :) = g_prev(num+1:end, :);
            g_est(end-num:end, :) = g_est_opt(curr_ind:curr_ind+num, :);
            out(:, 1+(j-1)*tau:j*tau) = g_est;
            g_prev = g_est;
            curr_ind = curr_ind+num+1;
        end
    end
end