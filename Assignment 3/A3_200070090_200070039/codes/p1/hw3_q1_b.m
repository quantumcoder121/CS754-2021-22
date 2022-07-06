True_Img = double(imread("barbara256.png")); %read image
Padded_Img = padarray(True_Img,[7, 7],0,"both"); %pad 0's to obtain smoother averaging
[m, n] = size(True_Img);
Reconst_Img_Padded = zeros(m+14,n+14);
Phi = normrnd(0,1,[32, 64]); %generate Phi (measurement matrix)

%we will use overlap of 7 pixels between consecutive patches of dim 8*8
%iterating over patches
for i = 1:263   
    for j = 1:263
        y=Phi*reshape(Padded_Img(i:i+7,j:j+7),[],1); %obtain measurement y
        x = reconst(y, Phi); 
        %add reconstructed patch to the appropriate position in the image
        Reconst_Img_Padded(i:i+7,j:j+7) = Reconst_Img_Padded(i:i+7,j:j+7)+reshape(x,8,8); 
    end
end
%unpad and average to obtain reconstructed image
Reconst_Img = Reconst_Img_Padded(8:263,8:263)/64; 
%rmse computed
rmse = norm(Reconst_Img-True_Img, 'fro')/norm(True_Img, 'fro');
imshow(uint8(Reconst_Img))

%This function reconstructs the vector x given y and Phi using ISTA
function x = reconst(y, Phi)
    U = kron(dctmtx(8)',dctmtx(8)');
    A = Phi*U;
    alpha = max(eig(A'*A))+1; %computing alpha using max eigenvalue
    theta = zeros(64,1); %initialise theta as 0
    for i = 1:200
        theta = wthresh(theta+A'*(y-A*theta)/alpha,'s',0.5/alpha); %update step
    end
    x = U*theta; %obtain x
end
