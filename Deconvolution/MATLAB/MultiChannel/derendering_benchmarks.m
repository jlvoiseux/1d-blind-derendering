num_lin_vec = linspace(10, 100, 10);
obs_size_move_vec = 2:4:40;
obs_size_source_vec = 2:4:40;
source_support_size_vec = linspace(1, 10, 10);
prop = 0.5;
tol = 1e-3;
alpha_vec = logspace(1e-1, 1e2, 10);
beta_vec = logspace(1e-3, 1, 10);

errx_num_lin_vec = zeros(1, 10);
errh_num_lin_vec = zeros(1, 10);
parfor i=1:10
    [errx_num_lin_vec(i), errh_num_lin_vec(i)] = DerenderingWrapper(round(num_lin_vec(i)), 20, 20, 4, prop, tol, 1, 0.01);    
    fprintf('Num Lin : Iteration %d finished.\n', i);
end
save("derendering_benchmark_num_lin.mat", 'errx_num_lin_vec', 'errh_num_lin_vec');

errx_obs_size_source = zeros(1, 10);
errh_obs_size_source = zeros(1, 10);
parfor i=1:10
    [errx_obs_size_source(i), errh_obs_size_source(i)] = DerenderingWrapper(40, 20, round(obs_size_source_vec(i)), 4, prop, tol, 1, 0.01);
    fprintf('Obs Size Source : Iteration %d finished.\n', i);
end
save("derendering_benchmark_obs_size_source.mat", 'errx_obs_size_source', 'errh_obs_size_source');

errx_alpha = zeros(1, 10);
errh_alpha = zeros(1, 10);
parfor i=1:10
    [errx_alpha(i), errh_alpha(i)] = DerenderingWrapper(40, 20, 20, 4, prop, tol, alpha_vec(i), 0.01);
    fprintf('Alpha : Iteration %d finished.\n', i);
end
save("derendering_benchmark_alpha.mat", 'errx_alpha', 'errh_alpha');

errx_beta = zeros(1, 10);
errh_beta = zeros(1, 10);
parfor i=1:10
    [errx_beta(i), errh_beta(i)] = DerenderingWrapper(40, 20, 20, 4, prop, tol, 1, beta_vec(i));
    fprintf('Beta : Iteration %d finished.\n', i);
end
save("derendering_benchmark_beta.mat", 'errx_beta', 'errh_beta');

errx_obs_size_move = zeros(1, 10);
errh_obs_size_move = zeros(1, 10);
parfor i=1:10
    [errx_obs_size_move(i), errh_obs_size_move(i)] = DerenderingWrapper(40, round(obs_size_move_vec(i)), 20, 4, prop, tol, 1, 0.01);
    fprintf('Obs Size Move : Iteration %d finished.\n', i);
end
save("derendering_benchmark_obs_size_move.mat", 'errx_obs_size_move', 'errh_obs_size_move');

errx_source_support_size = zeros(1, 10);
errh_source_support_size = zeros(1, 10);
parfor i=1:10
    [errx_source_support_size(i), errh_source_support_size(i)] = DerenderingWrapper(40, 20, 20, round(source_support_size_vec(i)), prop, tol, 1, 0.01);
    fprintf('Source Support Size : Iteration %d finished.\n', i);
end
save("derendering_benchmark_source_support_size.mat", 'errx_source_support_size', 'errh_source_support_size');




