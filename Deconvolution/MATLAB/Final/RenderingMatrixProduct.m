function out = RenderingMatrixProduct(s, brdf_mat, nmove, nsource, T)
   out = zeros(T, nsource, nmove);
    for j=1:nmove
        A = brdf_mat(:, :, j);    
        y = A*s;
        out(:, :, j) = y;       
    end
end

