epsilon_vals = [0.01, 0.1, 1, 10]; %values of epsilon
k_vals = [20, 25, 30, 35, 40, 45, 50, 55, 60]; %Sparsity values
prob = zeros(length(epsilon_vals)+1 ,length(k_vals)); %probability matrix
N_trials = 100; % number of trials
n_eps = length(epsilon_vals);
figure;
for k = 1: n_eps
    for j = 1: length(k_vals)
        for i = 1: N_trials
            [s_l1, s_rwl1] = rwl1_trial(epsilon_vals(k), k_vals(j));
            prob(k, j) = prob(k, j) + s_rwl1/N_trials;
            prob(n_eps+1, j) = prob(n_eps+1, j) + s_l1/(N_trials*n_eps);
        end
    end
    plot(k_vals, prob(k,:), '-o');
    hold on;
end
plot(k_vals, prob(n_eps+1,:), '-o');
hold on;
xlabel('Sparsity (k)');
ylabel('Recovery probability');
legend('epsilon = 0.01','epsilon = 0.1','epsilon = 1','epsilon = 10','unweighted L1');

function [succ_l1, succ_rwl1] = rwl1_trial(epsilon, n_non_zero)     
    % constants :
    n = 256; % size of x
    m = 100; % size of y
    n_iters = 5; % number of iterations of the algorithm

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
        w_diag_inv = abs(x_recon) + epsilon;
        
        %unweighted l1 reconstruction
        if j==1
            x_recon_l1 = x_recon;
        end
    end

    max_err_l1 = max(abs(x-x_recon_l1));
    if max_err_l1 < 1e-3
        succ_l1 = 1;
    else
        succ_l1 = 0;
    end

    max_err = max(abs(x-x_recon));
    if max_err < 1e-3
        succ_rwl1 = 1;
    else
        succ_rwl1 = 0;
    end
end