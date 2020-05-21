function [out, indices] = Sparsify(mat, prop)
    k = round(prop*numel(mat));
    [out, indices] = maxk(mat(:), k);
end

