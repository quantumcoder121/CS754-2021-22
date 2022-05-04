% constants :
n = 109; % size of x
N = n * n;
m = 3270; % size of y
n_iters = 5; % number of iterations of the algorithm
epsilon = 0.01; % correcting term for updating weights


% building x
img1 = imread("slice_50.png");
img2 = zeros(217, 217);
img2 = uint8(img2);
img2(19:199, 1:217) = img1;
img2 = double(img2);
img2 = imresize(img2, 0.5);
x = reshape(img2, N, 1);

% x = phantom('Shepp-Logan', 109);
% x = reshape(x, N, 1);

% Finding Radon Transform
% theta = 0:10:170;
% R = radon(img2, theta);

% Taking measurements y = Ax
Phi = randsrc(m, N);
Phi = orth(Phi')';
Phi_t = transpose(Phi);
y = Phi * x;


% reconstruction using reweighted l1 analysis

% initializing weights
w_diag = ones([N 1]);
% w_diag = ones([n 1]);

% running the algorithm
for j = 1:n_iters
    
    % w_diag_inv = zeros([n 1]);
    % for i = 1:n
        % w_diag_inv(i, 1) = 1 / (w_diag_inv(i, 1) + 0.01);
    % end
    % W = diag(w_diag);

    % reconstructing x with l1 magic package
    W = @(x) (w_diag .* dct_2d(x));
    Wt = @(x) (idct_2d(x) .* w_diag);
    % W = matrix_operator_2(N, w_diag);
    % Wt = W';
    opts.U = W;
    opts.Ut = Wt;
    x_recon = NESTA(Phi, Phi_t, y, 1e-3, 0.001, opts);

    % updating weights
    for i = 1:n
        w_diag(i, 1) = 1 / (abs(x_recon(i, 1)) + epsilon);
    end

end

% x_recon = reshape(uint8(normalize(x_recon, 'norm') * 255), 109, 109);
% img2_recon_2 = idct2(x_recon);
% img2_recon_2 = uint8(img2_recon_2);
% imshow(x_recon);

function ret = dct_2d(x)
n = 109;
x_mat = reshape(x, n, n);
x_dct = dct2(x_mat, n, n);
ret = reshape(x_dct, n * n, 1);
end

function ret = idct_2d(x)
n = 109;
x_mat = reshape(x, n, n);
x_idct = idct2(x_mat, n, n);
ret = reshape(x_idct, n * n, 1);
end