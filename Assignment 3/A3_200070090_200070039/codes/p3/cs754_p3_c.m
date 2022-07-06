% Reading and Converting Images
img1 = imread("slice_50.png");
img3 = imread("slice_51.png");
img2 = zeros(217, 217);
img4 = zeros(217, 217);
img2 = uint8(img2);
img4 = uint8(img4);
img2(19:199, 1:217) = img1;
img4(19:199, 1:217) = img3;
% imshow(img2);
% imshow(img4);
img2 = double(img2);
img4 = double(img4);

% Finding Radon Transform
theta = 0:10:170;
R = radon(img2, theta);
theta1 = 5:10:175;
R1 = radon(img4, theta);

% Reconstruction using Coupled Compressed Sensing
y = reshape(R, 1, []);
y1 = reshape(R1, 1, []);
y2 = horzcat(y, y1);
y2 = reshape(y2, [], 1);

m = 5562;
n = 47089;
A = coupled_matrix_operator(m, n, theta, theta1, 18);
lambda = 1;
rel_tol = 0.01;

[x, status] = l1_ls(A, y2, lambda, rel_tol);

x1 = x(1:n);
x2 = x(n + 1:end);
x2 = x2 + x1;
x1 = reshape(x1, 217, 217);
img2_recon_2 = idct2(x1);
img2_recon_2 = uint8(img2_recon_2);
imshow(img2_recon_2);
x2 = reshape(x2, 217, 217);
img4_recon_2 = idct2(x2);
img4_recon_2 = uint8(img4_recon_2);
% imshow(img4_recon_2);