k_vals = [20, 25, 30, 35, 40, 45, 50, 55, 60];
prob = zeros([length(k_vals),1]);
N_trials = 100;
figure;
for j = 1: length(k_vals)
    for i = 1: N_trials
        prob(j) = prob(j) + rwl1_adapt_eps_trial(k_vals(j))/N_trials;
    end
end
plot(k_vals, prob, '-o');
hold on;
xlabel('Sparsity (k)');
ylabel('Recovery probability');
legend('adaptive epsilon');

function succ = rwl1_adapt_eps_trial(n_non_zero)     
    % constants :
    n = 256; % size of x
    m = 100; % size of y
    n_iters = 5; % number of iterations of the algorithm
    i0 = int32(m/(4*log(n/m)));

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
        
        %adaptively assign epsilon
        x1 = maxk(x_recon,i0);
        epsilon = max(x1(i0), 1e-3);

        % updating weights    
        w_diag_inv = abs(x_recon) + epsilon;

    end

    max_err = max(abs(x-x_recon));
    if max_err < 1e-3
        succ = 1;
    else
        succ = 0;
    end
end