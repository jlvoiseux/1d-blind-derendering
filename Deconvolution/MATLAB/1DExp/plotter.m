% figure;
% subplot(1,2,1)
% plot(round(logspace(log10(50), log10(2000), 10)), errh_num_lin_vec, '-o');
% title("BRDF Estimation error");
% xlabel("N");
% ylabel("MSE");
% subplot(1,2,2)
% plot(round(logspace(log10(50), log10(2000), 10)), errx_num_lin_vec, '-o');
% title("Source Estimation error");
% xlabel("N");
% ylabel("MSE");
% 
% figure;
% subplot(1,2,1)
% plot(2:4:36, errh_obs_size_source(1:end-1), '-o');
% title("BRDF Estimation error");
% xlabel("B");
% ylabel("MSE");
% subplot(1,2,2)
% plot(2:4:36, errx_obs_size_source(1:end-1), '-o');
% title("Source Estimation error");
% xlabel("B");
% ylabel("MSE");
% 
% figure;
% subplot(1,2,1)
% plot(2:4:36, errh_obs_size_move(1:end-1), '-o');
% title("BRDF Estimation error");
% xlabel("A");
% ylabel("MSE");
% subplot(1,2,2)
% plot(2:4:36, errx_obs_size_move(1:end-1), '-o');
% title("Source Estimation error");
% xlabel("A");
% ylabel("MSE");
% 
% figure;
% subplot(1,2,1)
% imagesc(errh_obs_size_sourcemove);
% title("BRDF Estimation error");
% xlabel("A");
% ylabel("B");
% subplot(1,2,2)
% imagesc(errx_obs_size_sourcemove);
% title("Source Estimation error");
% xlabel("A");
% ylabel("B");

figure;
subplot(1,2,1)
imagesc(brdf_mat_flat);
title("True BRDF Matrix");
ylabel("N");
xlabel("4B");
subplot(1,2,2)
imagesc(brdf_est_flat);
title("Estimated BRDF Matrix");
ylabel("N");
xlabel("4B");

figure;
subplot(1,2,1)
imagesc(s);
title("True Source Matrix");
ylabel("M");
xlabel("B");
subplot(1,2,2)
imagesc(s_est);
title("Estimated Source Matrix");
ylabel("M");
xlabel("B");
