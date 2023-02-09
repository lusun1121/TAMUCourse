function [ elementToEdgeMap, nodeToEdgeMap] = extractElementsToEdges(n_nodes, elements, edges, minimalAngleInDegrees)
% return a list of element to edges with sign according to orientation
% of the global edge list.  That is,  if local edge ei is the same 
% orientation as the global edge index e, then we return the index e, 
% otherwise we return the negative index e.
%

% we first construct the node to edge map
nodeToEdgeMap = extractNodesToEdges(n_nodes,edges,minimalAngleInDegrees);

% Local orientation of element and edges
%
%  3
%  |\
%  | \
% e2  e1
%  |   \
%  |    \
%  1--e3-2

n_elements = length(elements(:,1));

% Construct the element to edge map
elementToEdgeMap = zeros(n_elements,3);
for i=1:n_elements
    edge1 = nonzeros(intersect(nodeToEdgeMap(elements(i,2),:),nodeToEdgeMap(elements(i,3),:)));
    edge2 = nonzeros(intersect(nodeToEdgeMap(elements(i,3),:),nodeToEdgeMap(elements(i,1),:)));
    edge3 = nonzeros(intersect(nodeToEdgeMap(elements(i,1),:),nodeToEdgeMap(elements(i,2),:)));
    
    %compare edgei to global ei orientation.  Positive if local edge is
    %same as global edge,  negative otherwise
    if (edges(edge1,1) == elements(i,2))
        elementToEdgeMap(i,1) = edge1;
    else
        elementToEdgeMap(i,1) = -edge1;
    end
    
    if (edges(edge2,1) == elements(i,3))
        elementToEdgeMap(i,2) = edge2;
    else
        elementToEdgeMap(i,2) = -edge2;
    end
    
    if ( edges(edge3,1) == elements(i,1))
        elementToEdgeMap(i,3) = edge3;
    else
        elementToEdgeMap(i,3) = -edge3;
    end

end

end