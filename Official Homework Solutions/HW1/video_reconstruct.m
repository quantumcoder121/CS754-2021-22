clear;
clc;

A = mmread('cars.avi');
T = 3;
for i=1:T
   X(:,:,i) = double(rgb2gray(A.frames(i).cdata)); 
end
[H,W,T] = size(X);
X = X(175:H,130:W,:);

C = rand(size(X));
C(C <= 0.5) = 0; C(C > 0.5) = 1;
Y = sum(X.*C,3);
Y = Y + randn(size(Y))*2;

v= Y; v(v < 0) = 0; imwrite(uint8(255*v/max(v(:))),sprintf('%d_coded_snapshot.jpg',T));

ps = 8; ps2 = ps*ps;
[H,W] = size(Y);
D = kron(dctmtx(ps)',dctmtx(ps)');

Xrec = zeros(size(X)); numcount = Xrec;
fprintf('\n');
for i=1:H-ps+1
    fprintf(' %d',i);
    for j=1:W-ps+1
        C_curr = C(i:i+ps-1,j:j+ps-1,:);
        Phi = zeros(ps2,ps2*T);
        for k=1:T
            c = squeeze(C_curr(:,:,k));
           Phi(:,(k-1)*ps2+1:k*ps2) = diag(c(:))*D; 
        end
        
        y = Y(i:i+ps-1,j:j+ps-1); y = y(:);
        thetas = omp_error(Phi,y,2); 
        
        for k=1:T
            Xrec(i:i+ps-1,j:j+ps-1,k) = Xrec(i:i+ps-1,j:j+ps-1,k) + reshape(D*thetas((k-1)*ps2+1:k*ps2),ps,ps);
            numcount(i:i+ps-1,j:j+ps-1,k) = numcount(i:i+ps-1,j:j+ps-1,k) + 1;
        end
    end
end
Xrec = Xrec./numcount;

for k=1:T
   imshow([Xrec(:,:,k) X(:,:,k)]/255); 
   
   v = [Xrec(:,:,k) X(:,:,k)]; 
   v(v < 0) = 0;
   v = uint8(v);
   fname = sprintf('%d_%d.png',T,k);
   imwrite(v,fname);
   
   pause; 
end

fprintf ('\nRMSE= %f',sum((X(:)-Xrec(:)).^2)/sum((X(:)).^2));
