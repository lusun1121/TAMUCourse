function [edgesToElementsMap, counts, nodesToElementsMap] = extractEdgesToElements(n_nodes, elements, edges, minimalAngleInDegrees) 
%
% The Edge to Element map.  each edge has one or two elements that it
% intersects.  For each edge, we return the index of the elements that 
% are it's neighbors and the count of how many neighbors and the by product
% nodeToElementMap can also be returned
%

% Local orientation of element and edges
%
%  3
%  |\
%  | \
% e2  e1
%  |   \
%  |    \
%  1--e3-2
%
% global edge has normal to the right which allows us to compute the jumps
% 
%   v2
%   |
%   |---> n
%   |
%   v1

nodesToElementsMap = extractNodesToElements(n_nodes,elements,minimalAngleInDegrees);

% sort edges smaller 
n_edges = size(edges,1);
edgesToElementsMap = zeros(n_edges,2);  % at most 2 elements per edge
counts = zeros(n_edges,1);

for i = 1:n_edges
    possibleElements = nonzeros(intersect(nodesToElementsMap(edges(i,1),:), nodesToElementsMap(edges(i,2),:)));
    counts(i) = size(possibleElements,1);
    
    % determine if local edge on each triangle is in same orientation as 
    % global edge so that we can determine normal to edge
    
    local_ele_order = elements(possibleElements(1),:);
    K = strfind([local_ele_order, local_ele_order(1)], edges(i,:));        
    if isempty(K)
        edgesToElementsMap(i,1) = -possibleElements(1);
        if counts(i)>1
            edgesToElementsMap(i,2) = possibleElements(2);
        end
    else
        edgesToElementsMap(i,1) = possibleElements(1);
        if counts(i)>1
            edgesToElementsMap(i,2) = -possibleElements(2);
        end
    end
end



