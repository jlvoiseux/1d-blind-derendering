function [f_est_final, h_est_final] = IterativeBlindDeconv(g, h_size, f_size, alpha, count)
    G = fft(g);
    G_size = length(G);
    %f_est = ones(1, f_size);
    f_est = g;
    h_est = zeros(1, h_size);
    h_est(ceil(size(h_est, 2)/2)) = 1;
    
    H_est = fft(h_est, G_size);
    for i=1:count
        F_est = fft(f_est, G_size);
        H_tilde = apply_fourier_constraints(G, F_est, H_est, alpha);
        h_tilde = ifft(H_tilde);
        h_est = abs(h_tilde(1:h_size));
        
        H_est = fft(h_est, G_size);
        F_tilde = apply_fourier_constraints(G, H_est, F_est, alpha);
        f_tilde = ifft(F_tilde);
        f_est = abs(f_tilde(1:f_size));
    end
    f_est_final = f_est;
    h_est_final = h_est;
    
end

function out = apply_fourier_constraints(A, B, C, alpha)
    out = (A.*conj(B))./(abs(B).^2+(alpha./(abs(C).^2)));
end
