function [T]=gettriangle(i)
if i==1
    T = TriangleReader('domainforp2.1.node','domainforp2.1.ele','domainforp2.1.edge');
elseif i==2
        T = TriangleReader('domainforp2.2.node','domainforp2.2.ele','domainforp2.2.edge');
else
           T = TriangleReader('domainforp2.3.node','domainforp2.3.ele','domainforp2.3.edge');
end
