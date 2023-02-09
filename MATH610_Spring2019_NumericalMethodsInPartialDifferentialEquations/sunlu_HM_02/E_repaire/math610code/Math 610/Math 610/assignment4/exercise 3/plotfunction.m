

for cell=1:T.n_elements
    dofIndices = DoFHandler.dofs(cell,:);
    %vertical [v1;v2;v3]
    vertices = T.nodes(T.elements(cell,:),:);

    X = vertices(:,1); 
    Y = vertices(:,2); 
    Z=zeros(3,1);
    for i=1:3
        Z(i)=uh(dofIndices(i));
    end

     patch(X,Y,Z,'w');
      view(3);
%      patch(X,Y,Z);

end