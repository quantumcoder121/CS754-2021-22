% Reading and Converting Images
img1 = imread("slice_50.png");
% img1 = imread("slice_51.png");
img2 = zeros(217, 217);
img2 = uint8(img2);
img2(19:199, 1:217) = img1;
% imshow(img2);
img2 = double(img2);

% Reconstruction using Compressed Sensing
% l1_ls

m = 5562;
n = 47089;
A = matrix_operator(m, n, theta, 18);
y = reshape(R, [], 1);
lambda = 0.01;
rel_tol = 0.01;

[x, status] = l1_ls(A, y, lambda, rel_tol);

x = reshape(x, 217, 217);
img2_recon_2 = idct2(x);
img2_recon_2 = uint8(img2_recon_2);
imshow(img2_recon_2);