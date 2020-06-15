\bm{o}clear all
close all

num_lin = 1000;
sigma = 1;
blurred_mirror_BRDF = @(angle_in, angle_out, sigma) (normpdf(angle_in+angle_out, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(angle_in+angle_out, 15, sigma)/normpdf(0, 0, sigma) + normpdf(angle_in+angle_out, -35, sigma)/normpdf(0, 0, sigma)  +  normpdf(angle_in+angle_out, 35, sigma)/normpdf(0, 0, sigma));
%margin = 15;
%mirror_BRDF = @(angle_in, angle_out, margin) (1*(abs(angle_in-angle_out) <= margin));

obs_pos = [0, -5];
obs_size_move = 20;
obs_size_source = 20;
obs_interval = [45, -45]; % Convention : Move from left to right
obs_width = 60; % Must be smaller than the interval
focus_dist = 2.5;
obs = BuildObs(obs_pos, obs_size_move, obs_size_source, obs_interval, obs_width, focus_dist);

source_pos = [0, -1];
source_support_width = 5;
source_support_size = 4;
source = BuildSource(source_pos, source_support_width, source_support_size, obs_size_source, false);
s = reshape(source(:, 3, :),[source_support_size, obs_size_source]);

[x_axis, brdf_mat] = CreateRenderingMatrixFromBRDFFinal(obs, source, blurred_mirror_BRDF, num_lin, sigma);
g = RenderingMatrixProduct(s, brdf_mat, obs_size_move, obs_size_source, num_lin);

% wall_points = ComputeWallPoints(obs, num_lin, obs_size_move);
% wall_points_ids = ComputeWallPointsIds(wall_points);
% test0 = Flatten(brdf_mat);
% test1 = BuildOptFromFlat(test0, num_lin, obs_size_move, source_support_size, wall_points);
% [test2, testind] = Sparsify(test1, 0.05);
% test3 = UnSparsify(test2, testind, size(test1, 1), size(test1, 2));
% test4 = BuildFlatFromOpt(test3, num_lin, obs_size_move, source_support_size, wall_points, wall_points_ids);

% figure;
% subplot(1, 2, 1)
% imagesc(test0);
% subplot(1, 2, 2)
% imagesc(test4);
% 
% figure;
% subplot(1, 2, 1)
% imagesc(test1);
% subplot(1, 2, 2)
% imagesc(test3);
% 
% figure;
% imagesc(test0-test4);
% disp(immse(test0, test4));

[s_est, brdf_est] = BlindDerendering(g, source_support_size, 1e-3, 10, obs, 'iter', brdf_mat, s, 0.05, 30);
brdf_mat_flat = Flatten(brdf_mat);
brdf_est_flat = Flatten(brdf_est);

