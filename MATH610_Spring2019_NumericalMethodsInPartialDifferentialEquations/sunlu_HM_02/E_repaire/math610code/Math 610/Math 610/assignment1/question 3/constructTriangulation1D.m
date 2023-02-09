
function [T] = constructTriangulation1D(L, num_elements)

T.n_elements = num_elements;
T.n_nodes = T.n_elements+1;
T.n_edges = T.n_nodes;
T.nodes = [linspace(0,pi/6,3/5*num_elements)';...
    linspace(pi/6+(pi/4-pi/6)/(1/5*num_elements-1),pi/4,1/5*num_elements)';...
    linspace(pi/4+(1-pi/4)/(1/5*num_elements),L,1/5*num_elements+1)'];
T.elements = [(1:(T.n_nodes-1))', (2:T.n_nodes)' ];
T.edges = (1:T.n_edges)'; 

T.edgeFlags = zeros(T.n_edges,1);
T.edgeFlags(1) = 1;
T.edgeFlags(end) = 1; 

