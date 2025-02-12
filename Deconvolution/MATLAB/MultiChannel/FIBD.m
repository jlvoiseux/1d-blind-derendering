function [s_interf_est, g_interf_est] = FIBD(d_interf, T, n, tau, tol, obs, empty_source, mirror_brdf, num_lin, margin)    
    % Initial guesses
      s_interf_est = zeros(2*T-1, tau);
      s_interf_est(T) = 1;
      g_interf_est = rand(2*tau-1, n*n);
%     for i = 1:n
%         for j = 1:n
%             index = get_col_num(i, j, n);
%             if i == j 
%                 g_interf_est(:, index) = 0;
%                 g_interf_est(tau, index) = rand();
%             end
%         end
%     end
    indices = get_indices(obs, empty_source, mirror_brdf, num_lin, margin, true, tau);
    alpha = [0 0];
    % Loop over decreasing alpha
    for i = 1:length(alpha)
        W1 = inf;
        W2 = inf;
        deltaW = inf;
        while deltaW > tol
            % 1. Update s
            % 1.1 Optimize W with s
            sa = optimvar('sa', 2*T-1, tau);
            Wfun = @(sa) computeW(sa, g_interf_est, d_interf, n, tau, 2*T-1, indices, alpha(i));
            Wexp = fcn2optimexpr(Wfun,sa);
            Wprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Wexp);
            Wprob.Constraints.cons1 = sa(T) == 1;            
            Wprob.Constraints.cons2 = sa(1:T-1, :) == flip(sa(T+1:end, :));
            Wprob.Constraints.cons3 = sa(:) >= 0; 
            W0.sa = s_interf_est;            
            [Wsol,Wfval,~,~] = solve(Wprob,W0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 1e5));
            % 1.2 Assignements
            s_interf_est = Wsol.sa;
            W1p = W1;
            W1 = Wfval;
            % 2. Update g
            % 2.1 Optimize W with g
            gij = optimvar('gij', 2*tau-1, n*n);
            Wfun = @(gij) computeW(s_interf_est, gij, d_interf, n, tau, 2*T-1, indices, alpha(i));
            Wexp = fcn2optimexpr(Wfun,gij);
            Wprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Wexp); 
            W0.gij = g_interf_est;            
            [Wsol,Wfval,~,~] = solve(Wprob,W0,'Options', optimoptions(@fminunc,'MaxFunctionEvaluations', 1e5));
            % 2.2 Assignements
            g_interf_est = Wsol.gij;
            W2p = W2;
            W2 = Wfval;
            deltaW = max([W1p - W1 W2p - W2]);            
            disp(deltaW);
        end
    end
    g_interf_est = test_sum(2*tau-1, d_interf);
end

function out = get_col_num(i,j, obs_num)
    out = (i-1).*obs_num + j;
end

function out = computeW(s_interf_est, g_interf_est, d_interf, n, tau, T, indices, alpha)
    % 1.1.1 Compute V
    V = 0;
    s_interf_est = s_interf_est(:, 1);
    for k=1:n
        for l=1:n
            temp_diff = d_interf(:, get_col_num(k, l, n)) - conv(s_interf_est, fit_impulses(T, g_interf_est(:, get_col_num(k, l, n)),indices), 'same');
            temp_diff = temp_diff .* temp_diff;
            V = V + sum(temp_diff);
        end
    end
    % 1.1.2 Compute next term
    g_same = g_interf_est(:, get_col_num(1:n, 1:n, n));
    t = (-tau+1:tau-1)';
    Wterm = sum(sum(t.*t.*g_same.*g_same));
    W = V + alpha*Wterm;    
    out = W;
end

function out = test_sum(tau, interf_g)
    interf_g  = interf_g(5:45, :);
    target = zeros(tau, length(interf_g(1, :)));
    for i=1:length(interf_g(1, :))
        for j=1:length(target(:, 1))
            target(j, i) = mean(interf_g((j-1)*floor(length(interf_g(:, 1))/length(target(:, 1)))+1:j*floor(length(interf_g(:, 1))/length(target(:, 1))),i));
        end
    end
    out = target;
end

function indices = get_indices(obs, empty_source, mirror_brdf, num_lin, margin, interf, tau)
    g_mirror = SimpleRendering(obs, empty_source, mirror_brdf, num_lin, margin);
    g_mirror = g_mirror(1, :);
    g_mirror = trim(g_mirror, 2*tau-1);
    if interf
        g_mirror = xcorr(g_mirror); 
        g_mirror = trim(g_mirror, 2*tau-1);
        indices = find(g_mirror>0.5);
    else
        indices = find(g_mirror>0.5);
    end
end

function out = trim(signal, num)
    [~, ind] = maxk(signal, num);
    out = zeros(size(signal));
    out(ind) = 1;
end


function out = fit_impulses(T, impulses, indices)
    out = zeros(T, 1);
    out(indices) = impulses;
end
