function [u]=u_exact(x)
n=length(x);
u=zeros(size(x));
for i=1:n
    if 0.5<x(i) && x(i)<1
        u(i)=1;
    end
end