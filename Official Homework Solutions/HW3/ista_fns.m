function [x,J] = ista_fns(y,H,Ht,lambda,alpha,Nit)
% [x, J] = ista_fns(y, H, lambda, alpha, Nit)
% L1-regularized signal restoration using the iterated
% soft-thresholding algorithm (ISTA)
% Minimizes J(x) = norm2(y-H*x)^2 + lambda*norm1(x)
% INPUT
% y - observed signal
% H - function handle
% Ht - function handle for H’
% lambda - regularization parameter
% alpha - need alpha >= max(eig(H’*H))
% Nit - number of iterations
% OUTPUT
% x - result of deconvolution
% J - objective function
J = zeros(1, Nit); % Objective function
x = 0*Ht(y); % Initialize x
T = lambda/(2*alpha);
for k = 1:Nit
Hx = H(x);
J(k) = sum(abs(Hx(:)-y(:)).^2) + lambda*sum(abs(x(:)));
x = soft(x + (Ht(y - Hx))/alpha, T);
end
