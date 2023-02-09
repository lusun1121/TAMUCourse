T = TriangleReader('LShapedDomain.1.node','LShapedDomain.1.ele','LShapedDomain.1.edge');
% u=@(x)abs(x(1)-x(2));
u=@(x)sin(pi*(x(1)+x(2)));
%u=@(x)double(x(1)^2+x(2)^2>=1/4);
%For every elements
for ele=1:T.n_elements
    %vertical [v1;v2;v3]
    v = T.nodes(T.elements(ele,:),:);
  
    X = v(:,1);
    Y = v(:,2);
    Z=zeros(3,1);
    for i=1:3
        Z(i)=u(v(i,:));
    end

     patch(X,Y,Z,'w');
      view(3);
%     patch(X,Y,Z);
%      colorbar
end