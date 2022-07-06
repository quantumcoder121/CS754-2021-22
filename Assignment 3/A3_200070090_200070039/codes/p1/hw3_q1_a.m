True_Img = double(imread("barbara256.png")); %read image
[m, n] = size(True_Img);
%generate corrupted image. Noise has variance 3
Corrupted_Img = True_Img + normrnd(0,sqrt(3),[m, n]);
Reconst_Img = zeros([m,n]);
%iterate over patches
for i = 1:8:249
    for j = 1:8:249
        %obtain reconstructed patches
        Reconst_Img(i:i+7,j:j+7)=reconst(Corrupted_Img(i:i+7,j:j+7));
    end
end
%compute rmse of corrupted and reconstructed images
rmse_corr = norm(Corrupted_Img-True_Img, 'fro')/norm(True_Img, 'fro');
rmse_rec = norm(Reconst_Img-True_Img, 'fro')/norm(True_Img, 'fro');
imshow(uint8(Reconst_Img))

%This function denoises a patch Y using ISTA and gives X as reconstruction
function X = reconst(Y)
    % as matrix A is unitary, max eigenvalue of A'A is 1
    alpha = 2;  % alpha chosen accordingly
    theta = zeros(8,8); %initialise theta as 0
    for i=1:50
        theta = wthresh(theta+dct2(Y-idct2(theta))/alpha,'s',0.5/alpha); %update step
    end
    X = idct2(theta); %obtain X
end