
clear all; close all; clc;

% run triangle with the following flags on the A.poly file
%
% triangle -pq25a0.01e A.poly
%
% which generates A.1.{node, ele, edge} files
%

% T = TriangleReader('A.1.node','A.1.ele');
 T = TriangleReader('A.1.node','A.1.ele','A.1.edge');
% T = TriangleReader('SRP.1.node','SRP.1.ele','SRP.1.edge');

%T = TriangleReader('LShapedDomain.1.node','LShapedDomain.1.ele','LShapedDomain.1.edge');
%T = TriangleReader('LShapedDomain.2.node','LShapedDomain.2.ele','LShapedDomain.2.edge');
%T = TriangleReader('LShapedDomain.3.node','LShapedDomain.3.ele','LShapedDomain.3.edge');
%T = TriangleReader('LShapedDomain.4.node','LShapedDomain.4.ele','LShapedDomain.4.edge');

% plot elements fill blue
for i = 1:T.n_elements
    X = T.nodes( T.elements(i,:),1);
    Y = T.nodes( T.elements(i,:),2);
    patch(X,Y,'b');
end

figure
% another approach to plotting elements with no fill
for i = 1:T.n_elements
    V = T.nodes( T.elements(i,:), :);
    f = [1 2 3];
    patch('Faces', f, 'Vertices',V,...
       'EdgeColor','green','FaceColor','none','LineWidth',1);
end

figure
% another approach to plotting elements with fill color specified at each node
for i = 1:T.n_elements
    V = T.nodes( T.elements(i,:), :);
    f = [1 2 3];
    Col = V(:,1); % color by x coordinate
    patch('Faces', f, 'Vertices',V,'FaceVertexCData', Col,...
       'EdgeColor','black','FaceColor','interp','LineWidth',1);
end

% plot edges
figure
for i = 1:T.n_edges
    X = T.nodes( T.edges(i,:),1);
    Y = T.nodes( T.edges(i,:),2);
    line(X,Y)
end

% plot boundary edges
figure
for i = 1:T.n_boundaryedges
    X = T.nodes( T.edges(T.boundaryedges(i),:),1);
    Y = T.nodes( T.edges(T.boundaryedges(i),:),2);
    line(X,Y)
end

% plot interior edges
figure
for i = 1:T.n_interioredges
    X = T.nodes( T.edges(T.interioredges(i),:),1);
    Y = T.nodes( T.edges(T.interioredges(i),:),2);
    line(X,Y)
end

%title('h=0.1');
%title('h=0.2');
%title('h=0.05');
%title('h=0.025');

