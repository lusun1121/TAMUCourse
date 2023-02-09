function [nodeToElementMap, counts] = extractNodesToElements(n_nodes,elements,minimalAngleInDegrees)
%node to edge
if nargin < 3
    minimalAngleInDegrees = 28; % default to 28 degrees
end

% max num of triangles that could touch node based on 28 degree minimum
% angle.  The number of edges is the same as the number of elements
maxelements = ceil(360/minimalAngleInDegrees); 
n_ele = length(elements(:,1));

% make room for all edges that touch a node
nodeToElementMap=zeros(n_nodes,maxelements); % pointers to elements
counts=ones(n_nodes,1); % count of elements
for i=1:n_ele
    nodeToElementMap(elements(i,1),counts(elements(i,1)))=i; counts(elements(i,1))=counts(elements(i,1))+1;
    nodeToElementMap(elements(i,2),counts(elements(i,2)))=i; counts(elements(i,2))=counts(elements(i,2))+1;
    nodeToElementMap(elements(i,3),counts(elements(i,3)))=i; counts(elements(i,3))=counts(elements(i,3))+1;
end