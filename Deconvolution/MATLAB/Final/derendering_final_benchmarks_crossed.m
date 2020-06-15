disp("Derendering v2.1");
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

errx_obs_size_sourcemove = zeros(iter_avg, 10, 10);
errh_obs_size_sourcemove = zeros(iter_avg, 10, 10);
parfor i=1:iter_avg
    fprintf('Thread %d started.\n', i);
    obs_size_source_vec = 2:2:20;
    obs_size_move_vec = 2:2:20;
    count = 0;
    for j=1:10
        for k=1:10
            count = count + 1;
            fprintf('Thread %d : Progress = %d \n', i, count);
            [errx_obs_size_sourcemove(i, j, k), errh_obs_size_sourcemove(i, j, k)] = DerenderingBench(num_lin_standard, angle_standard, sigma_standard, obs_size_move_vec(k), obs_size_source_vec(j), sss_standard, alpha_standard, prop_standard, max_count, tol);    
        end
    end
end
errx_obs_size_sourcemove = reshape(mean(errx_obs_size_sourcemove, 1), [10 10]);
errh_obs_size_sourcemove = reshape(mean(errh_obs_size_sourcemove, 1), [10 10]);
save("derendering_benchmark_obs_size_sourcemove.mat", 'errx_obs_size_sourcemove', 'errh_obs_size_sourcemove');






