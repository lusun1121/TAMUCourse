function [a,b,k,l,g]=getfandinterval(v1,v2,num1,num2)
if v1(1)==v2(1)
    if v1(2)>v2(2)
        a=v2(2);
        b=v1(2);
        k=num2;
        l=num1;
    else
        a=v1(2);
        b=v2(2);
        k=num1;
        l=num2;
    end
    if v1(1)<=0
        g=@(x)2*pi*sin(pi*x);
    else
        g=@(x)-2*pi*sin(pi*x);
    end
else
    if v1(1)<v2(1)
        a=v1(1);
        b=v2(1);
        k=num1;
        l=num2;
    else
        a=v2(1);
        b=v1(1);
        k=num2;
        l=num1;
    end
    
    if v1(2)>=0
        g=@(x)-2*pi*sin(pi*x);
    else
        g=@(x)2*pi*sin(pi*x);
    end
end