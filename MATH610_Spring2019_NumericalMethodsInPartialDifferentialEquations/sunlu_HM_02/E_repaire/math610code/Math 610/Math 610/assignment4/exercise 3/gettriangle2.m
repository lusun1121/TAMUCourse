function [T]=gettriangle2(i)
if i==1
    T = TriangleReader('LShapedDomain1.1.node','LShapedDomain1.1.ele','LShapedDomain1.1.edge');
elseif i==2
        T = TriangleReader('LShapedDomain1.2.node','LShapedDomain1.2.ele','LShapedDomain1.2.edge');
elseif i==3
            T = TriangleReader('LShapedDomain1.3.node','LShapedDomain1.3.ele','LShapedDomain1.3.edge');
        else
            T = TriangleReader('LShapedDomain1.4.node','LShapedDomain1.4.ele','LShapedDomain1.4.edge'); 
end
