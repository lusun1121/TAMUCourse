function [T]=gettriangle(i)
if i==1
    T = TriangleReader('domainforp1.1.node','domainforp1.1.ele','domainforp1.1.edge');
elseif i==2
        T = TriangleReader('domainforp1.2.node','domainforp1.2.ele','domainforp1.2.edge');
else
           T = TriangleReader('domainforp1.3.node','domainforp1.3.ele','domainforp1.3.edge');
end
