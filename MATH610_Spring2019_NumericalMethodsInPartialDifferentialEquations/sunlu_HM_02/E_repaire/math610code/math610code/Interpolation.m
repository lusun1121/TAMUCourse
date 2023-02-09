function T=Interpolation(choose)
T = TriangleReader(['LShapedDomain.',num2str(choose),'.node'],['LShapedDomain.',num2str(choose),'.ele'],['LShapedDomain.',num2str(choose),'.edge']);
% u=@(x)abs(x(1)-x(2));
u=@(x)sin(pi*(x(1)+x(2)));
%u=@(x)double(x(1)^2+x(2)^2>=1/4);
%For every elements
figure;
for ele=1:T.n_elements
    %vertical [v1;v2;v3]
    V = T.nodes(T.elements(ele,:),:);
  
    X = V(:,1);
    Y = V(:,2);
    Z=zeros(3,1);
    for i=1:3
        Z(i)=u(V(i,:));
    end
    patch('XData',X,'YData',Y,'ZData', Z,'FaceVertexCData', Z,...
       'EdgeColor','black','FaceColor','interp','LineWidth',1);% color by Z coordinate
    view(3);
end
end