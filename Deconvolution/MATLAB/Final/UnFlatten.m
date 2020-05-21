function out = UnFlatten(array, dim2, dim3)
    dim1 = size(array, 1);
    out = reshape(array, [dim1, dim2, dim3]);
end

