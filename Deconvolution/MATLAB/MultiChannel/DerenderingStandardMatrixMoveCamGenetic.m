function [s_est_full, g_est_deflat, g_est_flat] = DerenderingStandardMatrixMoveCamGenetic(d_full, s_est_full, g_est, T, nmove, nsource, tau, tol)
    g_est_flat = reshape(g_est, [T, tau*nmove]);
    d_full_flat = reshape(d_full, [T, nmove*nsource]);
    ac_mat = buildAutocorrMat(d_full, T, nsource, nmove);
    ac_mat = round(mean(ac_mat, 2));
    R1 = inf;
    R2 = inf;
    deltaR = inf;
    alpha = 1;
    g_est_opt = buildOptFromFlat(g_est_flat, T, nmove, tau, ac_mat);
    figure;
    imagesc(g_est_flat);
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
    while deltaR > tol
        % 1. Update g
        % 1.1 Optimize W with g         
        gaoptions = optimoptions('ga','UseParallel',true,'Display','iter', 'InitialPopulationMatrix', g_est_opt(:)');
        Rfun = @(gij) computeROptG(s_est_full(:), gij, d_full_flat, alpha, tau, nmove, nsource, ac_mat, T);
        [g_est_opt_line, fval] = ga(Rfun, size(g_est_opt, 1)*tau, [],[],[],[],[],[],[],gaoptions);
        g_est_opt = reshape(g_est_opt_line, [length(g_est_opt_line)/tau, tau]);
        % 2.2 Assignements
        R2p = R2;
        R2 = fval;   
        % 2. Update s
        % 2.1 Optimize W with s
        gaoptions = optimoptions('ga','UseParallel',true,'Display','iter', 'InitialPopulationMatrix', s_est_full(:)');
        Rfun = @(sa) computeROptG(sa, g_est_opt(:), d_full_flat, alpha, tau, nmove, nsource, ac_mat, T);
        [s_line, fval] = ga(Rfun, tau*nsource, [],[],[],[],[],[],[],gaoptions);
        s_est_full = reshape(s_line, [tau, nsource]);        
        R1p = R1;
        R1 = fval;              
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
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 1e6, 'Display', 'iter', 'UseParallel', true, 'Algorithm', 'sqp'));
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
        %s_est_full = Rsol.sa;
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

function out = computeROptG(s_line, g_line, d_full_flat, alpha, tau, nmove, nsource, ac_mat, T)
    % 1.1.1 Compute V
    R = 0;
    g_est_opt = reshape(g_line, [length(g_line)/tau, tau]);
    s_est_full = reshape(s_line, [tau, nsource]);
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
    out = zeros(T+sum(ac_mat)+length(ac_mat), tau);
    init = round(nmove/2);
    num = round(mean(ac_mat));
    for j=init:nmove
        g_est = g_est_flat(:, 1+(j-1)*tau:j*tau);
        if j==init
            out(1:T, :) = g_est;
            curr_ind = T+1;
        else
            %num = ac_mat(j-1);
            out(curr_ind:curr_ind+num-1, :) = g_est(end-num+1:end, :);
            curr_ind = curr_ind+num;
        end
    end
    for j=init-1:-1:1
        g_est = g_est_flat(:, 1+(j-1)*tau:j*tau);
        %num = ac_mat(j);
        out(curr_ind:curr_ind+num-1, :) = g_est(1:num, :);
        curr_ind = curr_ind+num;
    end
end

function out = buildFlatFromOpt(g_est_opt, T, nmove, tau, ac_mat)
    out = zeros(T, tau*nmove);
    init = round(nmove/2);
    num = round(mean(ac_mat));
    for j=init:nmove        
        if j==init
            out(:, 1+(j-1)*tau:j*tau) = g_est_opt(1:T, :);
            g_prev = g_est_opt(1:T, :);
            curr_ind = T+1;
        else
            %num = ac_mat(j-1);
            g_est = zeros(T, tau);
            g_est(1:end-num, :) = g_prev(num+1:end, :);
            g_est(end-num+1:end, :) = g_est_opt(curr_ind:curr_ind+num-1, :);
            out(:, 1+(j-1)*tau:j*tau) = g_est;
            g_prev = g_est;
            curr_ind = curr_ind+num;
        end
    end
    g_prev = g_est_opt(1:T, :);
    for j=init-1:-1:1
        %num = ac_mat(j);
        g_est = zeros(T, tau);
        g_est(num+1:end, :) = g_prev(1:end-num, :);
        g_est(1:num, :) = g_est_opt(curr_ind:curr_ind+num-1, :);
        out(:, 1+(j-1)*tau:j*tau) = g_est;
        g_prev = g_est;
        curr_ind = curr_ind+num;
    end
end