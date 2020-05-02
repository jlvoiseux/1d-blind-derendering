num_lin_vec = linspace(10, 100, 20);
obs_size_move_vec = linspace(2, 40, 20);
obs_size_source_vec = linspace(2, 40, 20);
source_support_size_vec = linspace(1, 20, 20);
prop = 0.5;
tol = 1e-3;
alpha_vec = linspace(1, 20, 20);

errx_num_lin_vec = zeros(1, 20);
errh_num_lin_vec = zeros(1, 20);
for i=1:20
    [errx_num_lin_vec(i), errh_num_lin_vec(i)] = DerenderingWrapper(round(num_lin_vec(i)), 20, 20, 4, prop, tol, 1);
end

save("derendering_benchmark_num_lin.mat", 'errx_num_lin_vec', 'errh_num_lin_vec');

errx_obs_size = zeros(20, 20);
errh_obs_size = zeros(20, 20);
for i=1:20
    for j=1:20
        [errx_obs_size(i,j), errh_obs_size(i,j)] = DerenderingWrapper(40, round(obs_size_move_vec(i)), round(obs_size_source_vec(j)), 4, prop, tol, 1);
    end
end

save("derendering_benchmark_obs_size.mat", 'errx_obs_size', 'errh_obs_size');

errx_source_support_size = zeros(1, 20);
errh_source_support_size = zeros(1, 20);
for i=1:20
    [errx_source_support_size(i), errh_source_support_size(i)] = DerenderingWrapper(40, 20, 20, round(source_support_size_vec(i)), prop, tol, 1);
end

save("derendering_benchmark_source_support_size.mat", 'errx_source_support_size', 'errh_source_support_size');

errx_alpha = zeros(1, 20);
errh_alpha = zeros(1, 20);
for i=1:20
    [errx_alpha(i), errh_alpha(i)] = DerenderingWrapper(40, 20, 20, 4, prop, tol, alpha_vec(i));
end

save("derendering_benchmark_alpha.mat", 'errx_alpha', 'errh_alpha');



