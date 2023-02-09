
numiteration=6;
n=20;
num=n*2.^(0:(numiteration-1));
c=zeros(numiteration,1);
for i=1:numiteration
    A=getmatrix(num(i));
%      A=1/num(i)*toeplitz([2,-1,zeros(1,num(i)-2)]);
    lambda1=abs(max((eig(A))));
    lambda2=abs(min((eig(A))));
    c(i)=lambda1/lambda2;
    fprintf('The condition number for matrix A is %f\n',c(i));
    if i>1
        rate=c(i)/c(i-1);
       fprintf('The condition number rate is %f\n',rate);
       h=1/n;
       rate/(h^2)
    end
end
