beta_vals = [0.01,0.1,1,10]; %set of epsilon values used
k_vals = [130, 150, 160, 170, 180, 200, 220, 250]; % number of errors
prob = zeros(length(beta_vals)+1 ,length(k_vals)); %probablity matrix
N_trials = 100; %number of trials
n = 128; % size of signal
m = 4*n; % size of encoded signal
n_eps = length(beta_vals);
figure;
for k = 1: n_eps
    for j = 1: length(k_vals)
        for i = 1: N_trials
            [s_l1, s_rwl1] = rwl1_trial(beta_vals(k), k_vals(j));
            prob(k, j) = prob(k, j) + s_rwl1/N_trials;
            prob(n_eps+1, j) = prob(n_eps+1, j) + s_l1/(N_trials*n_eps);
        end
    end
    plot(k_vals/m, prob(k,:), '-o');
    hold on;
end
plot(k_vals/m, prob(n_eps+1,:), '-o');
hold on;
xlabel('k/m');
ylabel('Recovery probability');
legend('beta = 0.01','beta = 0.1','beta = 1','beta = 10','unweighted');

function [succ_l1, succ_rwl1] = rwl1_trial(beta, err_spars)
    % constants :
    n = 128; % size of signal
    m = 4*n; % size of encoded signal
    n_iters = 5; % number of reweighted iterations
    
    x = randn([n 1]); %original signal to be sent
    % building matrix A and encoding codeword as Ax
    A = randn(m, n);
    A = normc(A);
    y = A * x; %uncorrupted y
    non_zero_set = randi([1 m], err_spars, 1);
    y(non_zero_set)=-y(non_zero_set); %corrupted codeword recieved
    epsilon = beta*std(y);
    
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