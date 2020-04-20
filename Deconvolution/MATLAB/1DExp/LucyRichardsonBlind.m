function [f_est, h_est] =  LucyRichardsonBlind(g, f_guess, h_guess, iter, iter_lucy)
    f_est = f_guess;
    h_est = h_guess;
    for i=1:iter
        h_est = deconvlucy(g, f_est, iter_lucy);
        f_est = deconvlucy(g, h_est, iter_lucy);
    end
end
