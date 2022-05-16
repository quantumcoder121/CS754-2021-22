m = 200;
n = 500;
l0_x = 18;
v_set_frac = 0.1;
lambda_vals = [0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1 2 5 10 15 20 30 50 100];
net_vals = 17;
x = zeros([n 1]);

s = rng;
non_zero_indices = randi([1 m], l0_x, 1);
rng(s);
for i = 1:l0_x
    non_zero_index = non_zero_indices(i, 1);
    s = rng;
    random_number = rand;
    rng(s);
    random_number = random_number * 1000.0;
    x(non_zero_index, 1) = random_number;
end
x = double(x);
s = rng;
phi = 2 * randi([0 1], m, n) - 1;
rng(s);
phi = double(phi);
phi = phi / sqrt(m);
phi_x = phi * x;
mod_phi_x = abs(phi_x);
sigma = 0.05 * mean(mod_phi_x);
s = rng;
eta = normrnd(0, sigma, m, 1);
rng(s);
y = phi_x + eta;
v_set_size = m * v_set_frac;
r_set_size = (1 - v_set_frac) * m;

lambda_vals = reshape(lambda_vals, [], 1);
ve_g_vals = zeros([net_vals 1]);
[l2_x, num_finite_x] = sumsqr(x);
l2_x = sqrt(l2_x);
rmse_vals = zeros([net_vals 1]);

for q = 1:net_vals
    s = rng;
    v_set_indices = randi([1 m], v_set_size, 1);
    rng(s);
    v_set = zeros([v_set_size 1]);
    v_set = double(v_set);
    for i = 1:v_set_size
        v_set_index = v_set_indices(i, 1);
        v_set(i, 1) = y(v_set_index, 1);
    end
    y_indices = 1:m;
    y_indices = reshape(y_indices, [], 1);
    r_set_indices = setdiff(y_indices, v_set_indices);
    r_set_indices = reshape(r_set_indices, [], 1);
    r_set = zeros([r_set_size 1]);
    for i = 1:r_set_size
        r_set_index = r_set_indices(i, 1);
        r_set(i, 1) = y(r_set_index, 1);
    end
    phi_r = zeros([r_set_size n]);
    for i = 1:r_set_size
        r_set_index = r_set_indices(i, 1);
        for j = 1:n
            phi_r(i, j) = phi(r_set_index, j);
        end
    end
    phi_v = zeros([v_set_size n]);
    for i = 1:v_set_size
        v_set_index = v_set_indices(i, 1);
        for j = 1:n
            phi_v(i, j) = phi(v_set_index, j);
        end
    end
    lambda = lambda_vals(q, 1);
    x_g = l1_ls(phi_r, r_set, lambda);
    v_diff = v_set - phi_v * x_g;
    [mse, num_finite] = sumsqr(v_diff);
    [l2_v_set, num_finite_2] = sumsqr(v_set);
    l2_v_set = sqrt(l2_v_set);
    ve_g = mse / l2_v_set;
    ve_g_vals(q, 1) = ve_g;
    x_diff = x_g - x;
    [l2_diff, num_finite_diff] = sumsqr(x_diff);
    l2_diff = sqrt(l2_diff);
    rmse = l2_diff / l2_x;
    rmse_vals(q, 1) = rmse;
end

lambda_vals = log(lambda_vals);

% plot(lambda_vals, ve_g_vals);
plot(lambda_vals, rmse_vals);