clear;
clc;

n = 500; 
m = 200;
s = 18;
x = zeros(n,1);
indices = randperm(n);
x(indices(1:s)) = rand(s,1)*1000;

Phi = rand(m,n); Phi(Phi <= 0.5) = -1; Phi(Phi > 0.5) = 1; Phi = Phi/sqrt(m);
y = Phi*x;
sigval = 0.05*mean(abs(y));
y = y + randn(m,1)*sigval;

Lambda = [0.0001 0.0005 0.001 0.005 0.01 0.05 0.1 0.5 1 2 5 8 9  10 11 12 15 20 30 50 100];
len = length(Lambda);

indices = randperm(m);
r_indices = indices(1:ceil(0.9*m));
v_indices = indices(ceil(0.9*m)+1:m);

ve = zeros(len,1);
for i=1:len
    fprintf('\n%f',Lambda(i));
    [x_est,status] = l1_ls(Phi(r_indices,:),y(r_indices),Lambda(i),0.01,true);
    ve(i) = sum((y(v_indices) - Phi(v_indices,:)*x_est).^2)/length(v_indices);
    rmseval(i) = norm(x_est-x,2)/norm(x,2);
    fprintf('\nstatus = %s, VE = %f, RMSE = %f',status,ve(i),rmseval(i));
end

plot(log(Lambda),ve/10000);
hold on; plot(log(Lambda),rmseval,'color','red');
xlabel('log(Lambda)');
legend('Validation error/10000','RMSE');

[minval,minind] = min(ve);
fprintf('\nMin. VE = %f, at Lambda = %f, RMSE at this value of Lambda = %f',minval,Lambda(minind(1)),rmseval(minind(1)));
[minval,minind] = min(rmseval);
fprintf('\nMin. RMSE = %f, at Lambda = %f',minval,Lambda(minind(1)));


