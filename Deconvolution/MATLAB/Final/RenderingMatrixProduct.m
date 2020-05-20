function out = RenderingMatrixProduct(s, brdf_mat, n, k, T)
   out = zeros(T, k, n);
    for j=1:n
        A = brdf_mat(:, :, j);    
        y = A*s;
        out(:, :, j) = y;       
    end
end

