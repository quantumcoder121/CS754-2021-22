% function for implementing transpose of matrix A
function x = ATMatrix(y, theta, sq_size, n_angles)
y = reshape(y, [], n_angles);
x_or = iradon(y, theta, 'linear', 'Ram-Lak', 1, sq_size);
x = dct2(x_or);
x = reshape(x, [], 1);
end