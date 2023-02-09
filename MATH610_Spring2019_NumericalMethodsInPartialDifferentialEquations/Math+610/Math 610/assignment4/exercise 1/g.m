function [g]=g(x)
if x(1)==-1 || x(1)==0
    g=@(x)2*pi*sin(pi*x(2));
elseif x(1)==-1
    g=@(x)-2*pi*sin(pi*x(2));
elseif x(2)==-1
    g=@(x)2*pi*sin(pi*x(1));
else
    g=@(x)-2*pi*sin(pi*x(1));
end