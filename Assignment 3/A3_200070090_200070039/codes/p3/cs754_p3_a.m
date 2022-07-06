% Reading and Converting Images
img1 = imread("slice_50.png");
% img1 = imread("slice_51.png");
img2 = zeros(217, 217);
img2 = uint8(img2);
img2(19:199, 1:217) = img1;
% imshow(img2);
img2 = double(img2);

% Finding Radon Transform
theta = 0:10:170;
R = radon(img2, theta);

% Reconstruction using Ram-Lak Filter
img2_recon_1 = iradon(R, theta, 'linear', 'Ram-Lak', 1, 217);
img2_recon_1 = uint8(img2_recon_1);
imshow(img2_recon_1);