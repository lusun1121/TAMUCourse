

 T = TriangleReader('LShapedDomain1.4.node','LShapedDomain1.4.ele','LShapedDomain1.4.edge');
for i = 1:T.n_elements
    X = T.nodes( T.elements(i,:),1);
    Y = T.nodes( T.elements(i,:),2);
    patch(X,Y,'w');
end
title('h = 0.025');
figure
% another approach to plotting elements with no fill
for i = 1:T.n_elements
    V = T.nodes( T.elements(i,:), :);
    f = [1 2 3];
    patch('Faces', f, 'Vertices',V,...
       'EdgeColor','green','FaceColor','none','LineWidth',1);
end

% figure
% % another approach to plotting elements with fill color specified at each node
% for i = 1:T.n_elements
%     V = T.nodes( T.elements(i,:), :);
%     f = [1 2 3];
%     Col = V(:,1); % color by x coordinate
%     patch('Faces', f, 'Vertices',V,'FaceVertexCData', Col,...
%        'EdgeColor','black','FaceColor','interp','LineWidth',1);
% end
% 
% % plot edges
% figure
% for i = 1:T.n_edges
%     X = T.nodes( T.edges(i,:),1);
%     Y = T.nodes( T.edges(i,:),2);
%     line(X,Y)
% end
% 
% % plot boundary edges
% figure
% for i = 1:T.n_boundaryedges
%     X = T.nodes( T.edges(T.boundaryedges(i),:),1);
%     Y = T.nodes( T.edges(T.boundaryedges(i),:),2);
%     line(X,Y)
% end
% 
% % plot interior edges
% figure
% for i = 1:T.n_interioredges
%     X = T.nodes( T.edges(T.interioredges(i),:),1);
%     Y = T.nodes( T.edges(T.interioredges(i),:),2);
%     line(X,Y)
% end
% 
