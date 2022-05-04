clear; clc;
v = VideoReader('cars.avi');
T = 3;
s = 105;
%reading frames and generating snapshot
frame= read(v,s);
frame = rgb2gray(frame);
[r1, c] = size(frame);
merged_frame = double(imcrop(frame, [c-300, r1-160, 240-1 ,120-1]));
C = double(randi([0,1],120, 240));
snapshot = C.*merged_frame;
C_mat = C;
i = s+1;
while i < T+s
    frame = read(v,i);
    frame = rgb2gray(frame);
    [r1, c] = size(frame);
    frame1 = double(imcrop(frame, [c-300, r1-160, 240-1 ,120-1]));
    C = double(randi([0,1],120, 240));
    snapshot = snapshot + C.*frame1; 
    C_mat = cat(2,C_mat,C);
    merged_frame = cat(1, merged_frame,frame1);
    i = i + 1;
end
%cropped frame dimensions
m = 120;
n =240;
sigma = 2; %std deviation of noise
snapshot = snapshot + normrnd(0,sigma,[m,n]); %corrupted by noise
p = 8; %patch dimension
ov=4; %patch overlap
N_vert = fix((m-ov)/(p-ov));
N_hor = fix((n-ov)/(p-ov));
Psi = kron (dctmtx(p)', dctmtx(p)');
l = 1;
k = 1;
x3=[];
%carrying out patchwise reconstruction
for p1=1:N_vert
    x2 = [];
    k=1;
    for p2=1:N_hor
        y = reshape(snapshot(l:(l+p-1),k:(k+p-1)),[],1);
        A = diag(reshape(C_mat(l:(l+p-1),k:(k+p-1)),[],1))*Psi;
        i = 1;
        while i < T
            A=cat(2, A, diag(reshape(C_mat(l:(l+p-1),(k+i*n):(k+p+i*n-1)),[],1))*Psi);
            i = i+1;
        end
        
        lambda = 50; %parameter for ISTA algo
        w_diag_inv = ones([p^2*T 1]); % initializing weights
        n_iters = 5; %number of iterations
        epsilon = 1;

        % running the algorithm
        for j = 1:n_iters
           
            % calculating matrices required for computation and initializing 'guess'
            W_inv = diag(w_diag_inv);
            A_dash = A * W_inv;
            init_guess = randn([p^2*T 1]); % our initial guess for x
            % reconstructing theta with ISTA Algo
            theta_dash = ISTA(A_dash, y , lambda);
            theta = W_inv * theta_dash;

            % updating weights    
            w_diag_inv = abs(theta) + epsilon;
        end
        theta = reshape(theta,p,p*T);
        x = idct2(theta(:,1:p));
        for i=2:T
            x = cat(3, x, idct2(theta(:,((i-1)*p+1):(i*p))));
        end
        if p2 == 1
            x2 = x;
        else
            x2(:,(end-ov+1):end,:)=(x2(:,(end-ov+1):end,:)+x(:,1:ov,:))/2;
            x2 = cat(2,x2,x(:,(end-p+ov+1):end,:));
        end
        k=k+p-ov;
    end
        if p1 == 1
            x3 = x2;
        else
            x3((end-ov+1):end,:,:)=(x3((end-ov+1):end,:,:)+x2(1:ov,:,:))/2;
            x3 = cat(1,x3,x2((end-p+ov+1):end,:,:));
        end
    l=l+p-ov;
end
recon_frames = uint8(x3);
%imshow(uint8(snapshot/T));
a1 = reshape(permute(recon_frames,[2,1,3]),size(recon_frames,2),[])';
a2 = uint8(merged_frame);
MSE= immse(a1,a2);
PSNR = 10*log10(256^2/MSE);
imshow(a2);
imshow(a1);

function theta = ISTA(A, y, lambda)
    alpha = max(eig(A'*A))+1; %computing alpha using max eigenvalue
    theta = zeros(size(A,2),1); %initialise theta as 0
    for i = 1:200
        theta = wthresh(theta+A'*(y-A*theta)/alpha,'s',0.5*lambda/alpha); %update step
    end
end
