function out = Flatten(array)
    [dim1, dim2, dim3] = size(array, [1 2 3]);
    out = reshape(array, [dim1, dim2*dim3]);
end

