function out = UnSparsify(vals, indices, dim1, dim2)
    out = zeros(dim1, dim2);
    out(indices) = vals;
end

