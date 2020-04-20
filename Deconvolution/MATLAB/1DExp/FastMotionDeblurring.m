function [f_est, h_est] = FastMotionDeblurring(g, f_size, h_size, iter, iter_shock, dt_shock, sigma1, sigma2, alpha, beta, thresh, thresh_fin)
    K = zeros(1, h_size);
    L = g;
    for i=1:iter
        L = L((end-f_size)/2+1:(end+f_size)/2);
        L = filter_bilateral(L, sigma1, sigma2);
        L = filter_shock(L, iter_shock, dt_shock);
        L = gradient_mag_thresholding(L, thresh);
        K = kernel_estimation(L, g, beta, h_size);
        L = deconvlucy(g, K, 10);
        L = L/norm(L);
        sigma2 = 0.9*sigma2;
        iter_shock = 0.9*iter_shock;
        thresh = min([1.1*thresh 1]);
    end
    f_est = abs(L((end-f_size)/2+1:(end+f_size)/2));
    f_est = f_est/norm(f_est);
    f_est = f_est/max(f_est);
    f_est(f_est < 0.75) = 0;
    h_est = abs(K);
    h_est = h_est/norm(h_est);
    h_est = h_est/max(h_est);
    h_est(h_est < 0.75) = 0;
end

function out = filter_bilateral(in, sigma1, sigma2)
    out = zeros(1, length(in));
    for p=1:length(in)
        gp = 0;
        w = 0;
        n = round(sqrt(sigma1) * 3);
        for i = -n:n
            q = max(1, min(length(in), p+i));
            g1 = normpdf(i, 0, sigma1);
            g2 = normpdf(in(p) - in(q), 0, sigma2);
            gp = gp + g1 * g2 * in(q);
            w = w + g1 * g2;
        end
        gp = gp ./ w;
        out(p) = gp;
    end
end

function out = filter_shock(in, iter, dt)
    out = in;
    for i=1:iter        
        out = out - sign(4*del2(in))*norm(gradient(in))*dt;
    end
end

function out = gradient_mag_thresholding(in, thresh)
    grad = abs(gradient(in));
    [max_e, max_i] = maxk(grad, round(thresh.*length(grad)));
    grad(setdiff(1:end,max_i)) = 0;
    out = grad;
end

function out = kernel_estimation(P, B, beta, h_size)    
    A = convmtx(P', h_size);
    b = B';
    out = (pcg(2*(A'*A) + 2*beta, 2*A'*b, 1e-6, 10000))';
    out = out./norm(out);
end

% clear all
% close all
% 
% x = linspace(-1, 1, 1000);
% mu = 0;
% sigma = 0.01;
% h_prep1 = normpdf(x, mu, sigma)/normpdf(0, mu, sigma);
% delta_num = 4;
% h_prep2 = zeros(1, 1000);
% h_prep2(1+length(h_prep2)/delta_num/2:length(h_prep2)/delta_num:end) = 1;
% h = conv(h_prep1, h_prep2);
% 
% 
% delta_num = 2;
% f = zeros(1, 1000);
% f(1+length(f)/delta_num/2:length(f)/delta_num:end) = 1;
% subplot(1, 2, 1)
% stem(f);
% title('Original signal')
% subplot(1, 2, 2);
% stem(h);
% title('Impulse response')
% 
% figure;
% g = conv(h, f);
% g = g/norm(g);
% stem(g);
% title('Convolved signal')
% 
% figure;
% [f_est, h_est] = FastMotionDeblurring(g, length(f), length(h), 10, 1, 1, 2, 0.5, 1, 5, 0.1, 0.1);
% subplot(1, 2, 1)
% stem(f_est);
% subplot(1, 2, 2);
% stem(h_est);