disp("Derendering v1.0");
num_workers = input("Number of workers (default : 8)");
parpool('local', num_workers);

iter_avg = input("Number of iterations for averaging (default : 16)");
max_count = input("Max Count (default : 30)");
tol = input("Tolerance (default : 1e-2)");

num_lin_standard = input("Num Lin Standard (default : 500)");
sigma_standard = input("Sigma Standard (default : 1)");
osm_standard = input("Obs Size Move Standard (default : 20)");
oss_standard = input("Obs Size Source Standard (default : 20)");
sss_standard = input("Source Support Size Standard (default : 4)");
alpha_standard = input("Alpha Standard (default : 0.1)");
prop_standard = input("Sparsity Proportion (default : 0.05)");
angle_standard = input("Standard Obs Angle (default : 60)");

errx_obs_size_source = zeros(iter_avg, 10);
errh_obs_size_source = zeros(iter_avg, 10);
parfor i=1:iter_avg
    fprintf('obs_size_source : Iteration %d started.\n', i);
    obs_size_source_vec = 2:4:40;
    for j=1:10
        fprintf('obs_size_source : Subiteration %d started.\n', j);
        [errx_obs_size_source(i, j), errh_obs_size_source(i, j)] = DerenderingBench(num_lin_standard, angle_standard, sigma_standard, osm_standard, obs_size_source_vec(j), sss_standard, alpha_standard, prop_standard, max_count, tol);     
    end
end
errx_obs_size_source = mean(errx_obs_size_source, 1);
errh_obs_size_source = mean(errh_obs_size_source, 1);
save("derendering_benchmark_obs_size_source.mat", 'errx_obs_size_source', 'errh_obs_size_source');

errx_alpha = zeros(iter_avg, 10);
errh_alpha = zeros(iter_avg, 10);
parfor i=1:iter_avg
    fprintf('alpha : Iteration %d started.\n', i);
    alpha_vec = logspace(log10(1e-2), log10(100), 10);
    for j=1:10
        fprintf('alpha : Subteration %d started.\n', j);
        [errx_alpha(i, j), errh_alpha(i, j)] = DerenderingBench(num_lin_standard, angle_standard, sigma_standard, osm_standard, oss_standard, sss_standard, alpha_vec(j), prop_standard, max_count, tol);     
    end
end
errx_alpha = mean(errx_alpha, 1);
errh_alpha = mean(errh_alpha, 1);
save("derendering_benchmark_alpha.mat", 'errx_alpha', 'errh_alpha');

errx_obs_size_move = zeros(iter_avg, 10);
errh_obs_size_move = zeros(iter_avg, 10);
parfor i=1:iter_avg
    fprintf('obs_size_move : Iteration %d started.\n', i);
    obs_size_move_vec = 2:4:40;
    for j=1:10
        fprintf('obs_size_move : Subteration %d started.\n', j);
        [errx_obs_size_move(i, j), errh_obs_size_move(i, j)] = DerenderingBench(num_lin_standard, angle_standard, sigma_standard, obs_size_move_vec(j), oss_standard, sss_standard, alpha_standard, prop_standard, max_count, tol);     
    end
end
errx_obs_size_move = mean(errx_obs_size_move, 1);
errh_obs_size_move = mean(errh_obs_size_move, 1);
save("derendering_benchmark_obs_size_move.mat", 'errx_obs_size_move', 'errh_obs_size_move');

errx_source_support_size = zeros(iter_avg, 10);
errh_source_support_size = zeros(iter_avg, 10);
parfor i=1:iter_avg
    fprintf('source_support_size : Iteration %d started.\n', i);
    source_support_size_vec = linspace(1, 10, 10);
    for j=1:10
        fprintf('source_support_size : Subiteration %d started.\n', j);
        [errx_source_support_size(i, j), errh_source_support_size(i, j)] = DerenderingBench(num_lin_standard, angle_standard, sigma_standard, osm_standard, oss_standard, source_support_size_vec(j), alpha_standard, prop_standard, max_count, tol);     
    end
end
errx_obs_size_move = mean(errx_obs_size_move, 1);
errh_obs_size_move = mean(errh_obs_size_move, 1);
save("derendering_benchmark_source_support_size.mat", 'errx_source_support_size', 'errh_source_support_size');

errx_sigma = zeros(iter_avg, 10);
errh_sigma = zeros(iter_avg, 10);
parfor i=1:iter_avg
    fprintf('sigma : Iteration %d started.\n', i);
    sigma_vec = linspace(0.1, 1, 10);
    for j=1:10
        fprintf('sigma : Subiteration %d started.\n', j);
        [errx_sigma(i, j), errh_sigma(i, j)] = DerenderingBench(num_lin_standard, angle_standard, sigma_vec(j), osm_standard, oss_standard, sss_standard, alpha_standard, prop_standard, max_count, tol);     
    end
end
errx_sigma = mean(errx_sigma, 1);
errh_sigma = mean(errh_sigma, 1);
save("derendering_benchmark_sigma.mat", 'errx_sigma', 'errh_sigma');

errx_prop = zeros(iter_avg, 10);
errh_prop = zeros(iter_avg, 10);
parfor i=1:iter_avg
    fprintf('prop : Iteration %d started.\n', i);
    prop_vec = linspace(0.01, 0.1, 10);
    for j=1:10
        fprintf('prop : Subiteration %d started.\n', j);
        [errx_prop(i, j), errh_prop(i, j)] = DerenderingBench(num_lin_standard, angle_standard, sigma_standard, osm_standard, oss_standard, sss_standard, alpha_standard, prop_vec(j), max_count, tol);     
    end
end
errx_prop = mean(errx_prop, 1);
errh_prop = mean(errh_prop, 1);
save("derendering_benchmark_prop.mat", 'errx_prop', 'errh_prop');

errx_angle = zeros(iter_avg, 10);
errh_angle = zeros(iter_avg, 10);
parfor i=1:iter_avg
    fprintf('angle : Iteration %d started.\n', i);
    angle_vec = linspace(15, 85, 10);
    for j=1:10
        fprintf('angle : Iteration %d started.\n', j);
        [errx_angle(i, j), errh_angle(i, j)] = DerenderingBench(num_lin_standard, angle_vec(j), sigma_standard, osm_standard, oss_standard, sss_standard, alpha_standard, prop_standard, max_count, tol);     
    end
end
errx_angle = mean(errx_angle, 1);
errh_angle = mean(errh_angle, 1);
save("derendering_benchmark_angle.mat", 'errx_angle', 'errh_angle');

errx_num_lin_vec = zeros(iter_avg, 10);
errh_num_lin_vec = zeros(iter_avg, 10);
parfor i=1:iter_avg
    fprintf('Num Lin : Iteration %d started.\n', i);
    num_lin_vec = round(logspace(log10(50), log10(2000), 10));
    for j=1:10
        fprintf('Num Lin : Subiteration %d started.\n', j);
        [errx_num_lin_vec(i, j), errh_num_lin_vec(i, j)] = DerenderingBench(num_lin_vec(j), angle_standard, sigma_standard, osm_standard, oss_standard, sss_standard, alpha_standard, prop_standard, max_count, tol);    
    end
end
errx_num_lin_vec = mean(errx_num_lin_vec, 1);
errh_num_lin_vec = mean(errh_num_lin_vec, 1);
save("derendering_benchmark_num_lin.mat", 'errx_num_lin_vec', 'errh_num_lin_vec');

