function [A]=omp_error(D,X,errorGoal) 

[~,P]=size(X);
[n,K]=size(D);
E2 = errorGoal^2*n; % error bound is on the residual squared norm
maxNumCoef = size(X,1);
A = zeros(size(D,2),size(X,2));

DN = D;
for i=1:K
    DN(:,i) = DN(:,i)/norm(DN(:,i)); 
end

for k=1:1:P
    if mod(k,500) == 0, fprintf ('%d ',k); end;

    x=X(:,k);
    residual=x;
	indx = zeros(1,n);
	a = [];
	
    currResNorm2 = sum(residual.^2);
	j = 0;
    
    while currResNorm2>E2 && j <= maxNumCoef,
		j = j+1;
         
        proj = DN'*residual;
        pos = find(abs(proj)==max(abs(proj)));
        pos = pos(1);
        indx(j) = pos;
        
        %a = pinv(D(:,indx(1:j)))*x;
        a = D(:,indx(1:j))\x;
        
        residual = x-D(:,indx(1:j))*a;
		currResNorm2 = sum(residual.^2);
   end;
   
   if (~isempty(indx))
       A(indx(1:j),k)=a;
   end
end;

