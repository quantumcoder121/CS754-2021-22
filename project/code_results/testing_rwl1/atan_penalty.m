% constants :
n = 256; % size of x
m = 100; % size of y
n_iters = 5; % number of iterations of the algorithm
epsilon = 10; % correcting term for updating weights
n_non_zero = 30; % measure of sparsity of x

% the following two sections need to be changed for images

% building x
x = zeros([n 1]);
non_zero_set = randi([1 n], n_non_zero, 1);
x(non_zero_set)=randn([n_non_zero 1]);


% building matrix A and taking measurements y = Ax
A = randn(m, n);
A = normc(A);
y = A * x;

% reconstruction using reweighted l1 minimization

% initializing weights
w_diag_inv = ones([n 1]);

% running the algorithm
for j = 1:n_iters
   
    % calculating matrices required for computation and initializing 'guess'
    W_inv = diag(w_diag_inv);
    A_dash = A * W_inv;
    init_guess = randn([n 1]); % our initial guess for x

    % reconstructing x with l1 magic package
    x_recon_dash = l1eq_pd(init_guess, A_dash, A_dash', y);
    x_recon = W_inv * x_recon_dash;
    
    % updating weights    
    w_diag_inv = x_recon.^2 + epsilon^2;

end
recon_err = norm(x-x_recon)/norm(x);
max_err = max(abs(x-x_recon));