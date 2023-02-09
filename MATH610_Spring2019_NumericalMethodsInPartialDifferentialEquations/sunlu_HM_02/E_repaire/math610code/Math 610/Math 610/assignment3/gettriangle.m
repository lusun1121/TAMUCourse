function [T]=gettriangle(i)
if i==1
    T = TriangleReader('LShapedDomain.1.node','LShapedDomain.1.ele','LShapedDomain.1.edge');
elseif i==2
        T = TriangleReader('LShapedDomain.2.node','LShapedDomain.2.ele','LShapedDomain.2.edge');
elseif i==3
            T = TriangleReader('LShapedDomain.3.node','LShapedDomain.3.ele','LShapedDomain.3.edge');
        else
            T = TriangleReader('LShapedDomain.4.node','LShapedDomain.4.ele','LShapedDomain.4.edge'); 
end
