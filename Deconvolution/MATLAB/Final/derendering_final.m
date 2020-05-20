clear all
%close all

num_lin = 100000;
sigma = 1;
blurred_mirror_BRDF = @(angle_in, angle_out, sigma) (normpdf(angle_in+angle_out, -15, sigma)/normpdf(0, 0, sigma)  +  normpdf(angle_in+angle_out, 15, sigma)/normpdf(0, 0, sigma) + normpdf(angle_in+angle_out, -35, sigma)/normpdf(0, 0, sigma)  +  normpdf(angle_in+angle_out, 35, sigma)/normpdf(0, 0, sigma));
margin = 15;
mirror_BRDF = @(angle_in, angle_out, margin) (1*(abs(angle_in-angle_out) <= margin));

obs_pos = [0, -5];
obs_size_move = 1;
obs_size_source = 10;
obs_interval = [45, -45]; % Convention : Move from left to right
obs_width = 90; % Must be smaller than the interval
focus_dist = 2.5;
obs = BuildObs(obs_pos, obs_size_move, obs_size_source, obs_interval, obs_width, focus_dist);

source_pos = [0, -0.1];
source_support_width = 5;
source_support_size = 2;
source = BuildSource(source_pos, source_support_width, source_support_size, obs_size_source, false);
s = reshape(source(:, 3, :),[source_support_size, obs_size_source]);

[x_axis, brdf_mat] = CreateRenderingMatrixFromBRDFFinal(obs, source, mirror_BRDF, num_lin, margin);
g = RenderingMatrixProduct(s, brdf_mat, obs_size_move, obs_size_source, num_lin);

figure;
for j=1:obs_size_move
    subplot(1, obs_size_move, j)
    stem(x_axis(j,:), g(:, 1, j));
end




