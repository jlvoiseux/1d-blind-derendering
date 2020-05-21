function [s_est, brdf_est] = DerenderingOpt(g, s_est, brdf_est, T, nmove, nsource, tau, tol, alpha, wall_points, wall_points_ids, doDisplay, propSparse)

    brdf_est_flat = Flatten(brdf_est);
    g_flat = Flatten(g);
    R1 = inf;
    R2 = inf;
    deltaR = inf;
    brdf_est_opt = BuildOptFromFlat(brdf_est_flat, T, nmove, tau, wall_points);
    [brdf_est_sparse, brdf_est_sparse_indices] = Sparsify(brdf_est_opt, propSparse);
    brdf_est_dim1 = size(brdf_est_opt, 1);
    brdf_est_dim2 = size(brdf_est_opt, 2);
    
    while deltaR > tol        
        % 1. Update s
        % 1.1 Optimize W with s
        sa = optimvar('sa', tau, nsource);
        Rfun = @(sa) ComputeROpt(sa, brdf_est_sparse, brdf_est_sparse_indices, brdf_est_dim1, brdf_est_dim2, g_flat, alpha, tau, nmove, nsource, wall_points, wall_points_ids, T);
        Rexp = fcn2optimexpr(Rfun,sa);
        Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp);
        Rprob.Constraints.cons1 = sa(:) >= 0;
        Rprob.Constraints.cons2 = sa(:) <= 1;
        R0.sa = s_est;            
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 1e5, 'Display', doDisplay, 'Algorithm', 'sqp'));
        % 1.2 Assignements
        s_est = Rsol.sa;
        R1p = R1;
        R1 = Rfval; 
        % 2. Update g
        % 2.1 Optimize W with g         
        brdfij = optimvar('brdfij', size(brdf_est_sparse, 1));
        Rfun = @(brdfij) ComputeROpt(s_est, brdfij, brdf_est_sparse_indices, brdf_est_dim1, brdf_est_dim2, g_flat, alpha, tau, nmove, nsource, wall_points, wall_points_ids, T);
        Rexp = fcn2optimexpr(Rfun, brdfij);
        Rprob = optimproblem('ObjectiveSense', 'minimize', 'Objective', Rexp); 
        Rprob.Constraints.cons1 = brdfij(:) >= 0;
        Rprob.Constraints.cons2 = brdfij(:) <= 1;
        R0.brdfij = brdf_est_sparse;            
        [Rsol,Rfval,~,~] = solve(Rprob,R0,'Options', optimoptions(@fmincon,'MaxFunctionEvaluations', 2e5, 'Display', doDisplay, 'Algorithm', 'sqp'));
        % 2.2 Assignements
        brdf_est_sparse = Rsol.brdfij;
        R2p = R2;
        R2 = Rfval;                     
        deltaR = max([R1p - R1 R2p - R2]);            
        disp(deltaR);
    end
    brdf_est_opt = UnSparsify(brdf_est_sparse, brdf_est_sparse_indices, brdf_est_dim1, brdf_est_dim2);
    brdf_est_flat = BuildFlatFromOpt(brdf_est_opt, T, nmove, tau, wall_points, wall_points_ids);
    brdf_est = UnFlatten(brdf_est_flat, tau, nmove);

end