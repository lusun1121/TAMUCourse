function [T]=gettriangle(i)
if i==1
    T = TriangleReader('domainforp3.1.node','domainforp3.1.ele','domainforp3.1.edge');
elseif i==2
        T = TriangleReader('domainforp3.2.node','domainforp3.2.ele','domainforp3.2.edge');
else
           T = TriangleReader('domainforp3.3.node','domainforp3.3.ele','domainforp3.3.edge');
end
