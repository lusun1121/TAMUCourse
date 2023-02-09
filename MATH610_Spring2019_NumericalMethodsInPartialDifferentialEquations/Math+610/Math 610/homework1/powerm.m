function [lambda]=powerm(M)
n=size(M);
n=n(1);
maxiter=10000;
b0=ones(n,1);
lambda0=norm(b0);
b0=b0/lambda0;
for i=1:maxiter
   b1=M*b0;
   lambda1=norm(b1);
   b1=b1/lambda1;
   if norm(lambda0-lambda1)<0.00001
       break
   end
   b0=b1;
   lambda0=lambda1;
end
lambda=lambda0;
end

