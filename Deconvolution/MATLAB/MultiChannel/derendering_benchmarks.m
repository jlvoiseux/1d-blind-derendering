num_lin_vec = linspace(10, 100, 10);
obs_size_move_vec = linspace(1, 20, 10);
obs_size_source_vec = linspace(1, 20, 10);
source_support_size_vec = linspace(1, 10, 10);
prop = 0.5;
tol = 1e-3;
alpha_vec = linspace(1, 20, 10);

errx_num_lin_vec = zeros(1, 10);
errh_num_lin_vec = zeros(1, 10);
for i=1:10
    [errx_num_lin_vec(i), errh_num_lin_vec(i)] = DerenderingWrapper(round(num_lin_vec(i)), 16, 16, 4, prop, tol, 1);
    save("derendering_benchmark_num_lin.mat", 'errx_num_lin_vec', 'errh_num_lin_vec');
end

errx_obs_size = zeros(10, 10);
errh_obs_size = zeros(10, 10);
for i=1:10
    for j=1:10
        [errx_obs_size(i,j), errh_obs_size(i,j)] = DerenderingWrapper(25, round(obs_size_move_vec(i)), round(obs_size_source_vec(j)), 4, prop, tol, 1);
        save("derendering_benchmark_obs_size.mat", 'errx_obs_size', 'errh_obs_size');
    end
end

errx_source_support_size = zeros(1, 10);
errh_source_support_size = zeros(1, 10);
for i=1:10
    [errx_source_support_size(i), errh_source_support_size(i)] = DerenderingWrapper(25, 16, 16, round(source_support_size_vec(i)), prop, tol, 1);
    save("derendering_benchmark_source_support_size.mat", 'errx_source_support_size', 'errh_source_support_size');
end

errx_alpha = zeros(1, 10);
errh_alpha = zeros(1, 10);
for i=1:10
    [errx_alpha(i), errh_alpha(i)] = DerenderingWrapper(25, 16, 16, 4, prop, tol, alpha_vec(i));
    save("derendering_benchmark_alpha.mat", 'errx_alpha', 'errh_alpha');
end




