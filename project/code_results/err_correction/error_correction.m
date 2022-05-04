% constants :
n = 128; % size of signal
m = 4*n; % size of encoded signal
n_iters = 5; % number of reweighted iterations
epsilon = 0.1; % correcting term for updating weights
err_spars = 100; % measure of sparsity of error

x = randn([n 1]); %original signal to be sent
% building matrix A and encoding codeword as Ax
A = randn(m, n);
A = normc(A);
y = A * x; %uncorrupted y
non_zero_set = randi([1 m], err_spars, 1);
y(non_zero_set)=-y(non_zero_set); %corrupted codeword recieved


% reconstruction using reweighted l1decode

% initializing weights
w_diag_inv = ones([n 1]);

% running the algorithm
for j = 1:n_iters
    % calculating matrices required for computation and initializing 'guess'
    W_inv = diag(w_diag_inv);
    A_dash = A * W_inv;
    init_guess = randn([n 1]); % our initial guess for x

    % reconstructing x with l1 magic package
    x_recon_dash = l1decode_pd(init_guess, A_dash, A_dash', y);
    x_recon = W_inv * x_recon_dash;
    
    % updating weights    
    w_diag_inv = abs(x_recon) + epsilon;
end
recon_err = norm(x-x_recon)/norm(x);
max_err = max(abs(x-x_recon));