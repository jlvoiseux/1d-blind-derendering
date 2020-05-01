function [s_est_full, g_est_deflat] = DerenderingStandardMatrixMoveCamLS(d_full, s_est_full, g_est, T, nmove, nsource, tau, tol)
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
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@lsqlin,'Display', 'iter'));
        % 1.2 Assignements
        s_est_full = Rsol.sa;
        R1p = R1;
        R1 = Rfval;
        % 2. Update g
        % 2.1 Optimize W with g
        gij = optimvar('gij', T, tau*nmove);
        Rfun = @(gij) computeR(s_est_full, gij, d_full_flat, alpha, tau, nmove, nsource, ac_mat);
        Rexp = fcn2optimexpr(Rfun,gij);
        Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp); 
        Rprob.Constraints.cons1 = gij(:) >= 0;
        Rprob.Constraints.cons2 = gij(:) <= 1;
        R0.gij = g_est_flat;            
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@lsqlin,'Display', 'iter'));
        % 2.2 Assignements
        g_est_flat = Rsol.gij;
        R2p = R2;
        R2 = Rfval;         
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
            if j~=1
                for k=1:tau
                    num = max(ac_mat(j-1, i), 1);
                    move_vec = zeros(num, 1);
                    %conv_vec = [g_est_prev(1:end-num, k); move_vec];
                    %R = R + alpha*sum(g_est(:, k)-conv_vec);
                    conv_vec = conv(g_est_prev(:, k), move_vec);
                    R = R + sum(g_est(:, k)-conv_vec(1:end-num+1));
                end
            end
            g_est_prev = g_est;
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
            out(j-1, i) = ind-T;
        end
    end
end