N_trials = 100;
recon_err = zeros([N_trials 1]);
recon_err_l1 = zeros([N_trials 1]);
for trial = 1: N_trials   
    % constants :
    n = 256; % size of x
    m = 100; % size of y
    n_iters = 5; % number of iterations of the algorithm
    epsilon = 1; % correcting term for updating weights
    n_non_zero = 10; % measure of sparsity of x
    sigma = 0.04;
    delta = sigma*(m+2^1.5*m^0.5)^0.5;
    
    % the following two sections need to be changed for images
    
    % building x
    x = zeros([n 1]);
    non_zero_set = randi([1 n], n_non_zero, 1);
    x(non_zero_set)=randn([n_non_zero 1]);
    
    
    % building matrix A and taking measurements y = Ax
    A = randn(m, n);
    A = normc(A);
    noise = normrnd(0,sigma,[m,1]);
    y = A * x + noise;
    
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
        x_recon_dash = l1qc_logbarrier(init_guess, A_dash, A_dash', y, delta);
        x_recon = W_inv * x_recon_dash;
        
        % updating weights    
        w_diag_inv = abs(x_recon) + epsilon;
        %unweighted reconstruction
        if j==1
            x_recon_l1 = x_recon;
        end
    end
    recon_err(trial) = norm(x-x_recon)/norm(x);
    recon_err_l1(trial) = max(abs(x-x_recon_l1));
end
histogram(recon_err);
%histogram(recon_err_l1);