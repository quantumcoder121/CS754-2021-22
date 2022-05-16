n = 128;
U = zeros(n);
%generation of random U
column_i = randn(n, 1); 
U(:, 1) = column_i ./ norm(column_i); 
for i = 2:n
    column_norm = 0;
    while column_norm < 1e-6
        column_i = randn(n, 1);
        column_i = column_i - U(:, 1:i - 1) * (U(:, 1:i - 1).' * column_i);
        column_norm = norm(column_i);
    end
    U(:, i) = column_i ./ column_norm;
end
mu = zeros(n,1);
m_vals = [40, 50, 64, 80, 100, 120];
alpha_vals = [0, 3];
figure;
for alpha=1:length(alpha_vals)
    lambda_vals = zeros(n,1);
    for i=1:128
        lambda_vals(i) = i^(-alpha_vals(alpha));  %diagonal values of covariance matrix
    end
    Sigma = U * diag(lambda_vals) * U'; %covariance matrix
    rmse_vals = zeros(length(m_vals), length(alpha_vals));
    rng(1);
    for i = 1:length(m_vals)
        for j = 1:15 %run a few times
            x = mvnrnd(mu ,Sigma ,1)'; %x generated according to required distribution
            phi = sqrt(1/m_vals(i)) * randn(m_vals(i), n); %phi genrated randomly
            stddev = 0.01 * mean(abs(phi * x));
            y = phi * x + stddev * randn(m_vals(i), 1);  % measurements y
            x_recon = (inv(phi'*phi + stddev^2 * (U * diag(1./lambda_vals) * U'))) * phi'*y; %reconstruction step
            rmse_vals(i, alpha) = rmse_vals(i, alpha) + norm(x_recon-x)/sqrt(128)/15;
        end
    end
    plot(m_vals, rmse_vals(:, alpha), '-o');
    hold on;
end
xlabel('no. of measurements (m)');
ylabel('RMSE');
legend('alpha = 0', 'alpha = 3');