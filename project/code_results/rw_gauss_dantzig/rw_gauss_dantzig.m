N_trials = 100;
recon_err = zeros([N_trials 1]);
recon_err_l1 = zeros([N_trials 1]);
n_corr = zeros([N_trials 1]);
n_corr_l1 = zeros([N_trials 1]);
for trial = 1: N_trials    
    % constants :
    n = 256; % size of x
    m = 100; % size of y
    n_iters = 10; % number of iterations of the algorithm
    epsilon = 0.1; % correcting term for updating weights
    n_non_zero = 10; % measure of sparsity of x
    sigma = 0.1;
    delta = 1;
    alpha = 0.25;
    
    % the following two sections need to be changed for images
    
    % building x
    x = zeros([n 1]);
    non_zero_set = randi([1 n], n_non_zero, 1);
    x(non_zero_set) = randn([n_non_zero 1]);
    x = x + sign(x);
    
    
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
    
        % dantzig step
        x_recon_dash = l1dantzig_pd(init_guess, A_dash, A_dash', y, delta);
        x_recon = W_inv * x_recon_dash;
        
        %identify support and refine estimate
        support = [];
        for i = 1:n 
            if abs(x_recon(i)) > alpha*sigma
                support(end+1) = i;
            end
        end
        x_recon = zeros([n 1]);
        A_sup = A(:,support);
        x_recon(support) = (A_sup'*A_sup)\(A_sup'*y);
        
        % updating weights    
        w_diag_inv = abs(x_recon) + epsilon;
         %unweighted gauss dantzig reconstruction
        if j==1
            x_recon_l1 = x_recon;
            n_corr_l1(trial) = size(intersect(support,non_zero_set),1);
        end
    end
    n_corr(trial) = size(intersect(support,non_zero_set),1);
    recon_err(trial) = norm(x-x_recon)/norm(x);
    recon_err_l1(trial) = max(abs(x-x_recon_l1));
end
histogram(recon_err);
histogram(recon_err_l1);