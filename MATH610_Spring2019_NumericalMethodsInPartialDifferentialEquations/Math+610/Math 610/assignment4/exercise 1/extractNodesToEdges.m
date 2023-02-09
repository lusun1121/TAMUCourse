function [ nodeToEdgeMap, counts ] = extractNodesToEdges(n_nodes,edges,minimalAngleInDegrees)
%node to edge
if nargin < 3
    minimalAngleInDegrees = 28; % default to 28 degrees
end

% max num of triangles that could touch node based on 28 degree minimum
% angle.  The number of edges is the same as the number of elements
maxelements = ceil(360/minimalAngleInDegrees); 
n_edge = length(edges(:,1));

% make room for all edges that touch a node
nodeToEdgeMap=zeros(n_nodes,maxelements); % pointers to edges
counts=ones(n_nodes,1); % count of edges
for i=1:n_edge
    nodeToEdgeMap(edges(i,1),counts(edges(i,1)))=i; counts(edges(i,1))=counts(edges(i,1))+1;
    nodeToEdgeMap(edges(i,2),counts(edges(i,2)))=i; counts(edges(i,2))=counts(edges(i,2))+1; 
end

