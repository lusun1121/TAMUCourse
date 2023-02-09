function [gradu]=grad_u_exact(x)
n=length(x);
gradu=zeros(size(x));
for i=1:n
if 0<=x(i) && x(i)<pi/6
    gradu(i)=12/pi;
elseif pi/6<=x(i) && x(i)<pi/4
    gradu(i)=6/pi;
elseif pi/4<=x(i) && x(i)<=1
    gradu(i)=4/pi;
end
end