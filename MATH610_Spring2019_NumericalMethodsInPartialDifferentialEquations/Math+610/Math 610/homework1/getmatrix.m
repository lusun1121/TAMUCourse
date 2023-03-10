function [A]=getmatrix(n)
A1=zeros(2*n+1,2*n+1);
A2=zeros(2*n+1,2*n+1);
localA1=[4,-3,-1;-3,36,-3;-1,-3,4];
localA2=[4,-6,2;-6,12,-6;2,-6,4];
for i=1:n
    A1((2*i-1:2*i+1),(2*i-1:2*i+1))= A1((2*i-1:2*i+1),(2*i-1:2*i+1))+localA1;
    A2((2*i-1:2*i+1),(2*i-1:2*i+1))= A2((2*i-1:2*i+1),(2*i-1:2*i+1))+localA2;
end
h=1/n;
A=h/30*A1+1/h*A2;
