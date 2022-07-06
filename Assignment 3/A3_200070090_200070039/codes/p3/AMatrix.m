% function for implementing matrix A
function y = AMatrix(x, theta, sq_size)
x = reshape(x, sq_size, sq_size);
x_it = idct2(x);
R_func = radon(x_it, theta);
y = R_func;
y = reshape(y, [], 1);
end