function estimate = TikDeconv(g, h, alpha)
    H = convmtx(h, length(g));
    temp1 = H'*H + alpha*speye(size(H, 2));
    temp2 = H'*g';
    estimate = temp1\temp2;
end

