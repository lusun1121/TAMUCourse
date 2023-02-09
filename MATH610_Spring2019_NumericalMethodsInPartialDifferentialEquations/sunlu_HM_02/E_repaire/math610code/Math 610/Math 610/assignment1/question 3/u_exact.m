function [u]=u_exact(x)
n=length(x);
u=zeros(size(x));
for i=1:n
if 0<=x(i) && x(i)<=pi/6
    u(i)=12/pi*x(i);
elseif pi/6<x(i)&& x(i)<=pi/4
    u(i)=6/pi*x(i)+1;
elseif pi/4<x(i) && x(i)<=1
    u(i)=4/pi*x(i)+3/2;
end
end
