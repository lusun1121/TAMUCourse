%  Written by Spencer Patty
%  Jan 20, 2017
%  for use in Math 610 Lab
function [T] = constructTriangulation1D(L, num_elements)

T.n_elements = num_elements;
T.n_nodes = T.n_elements+1;
T.n_edges = T.n_nodes;
T.nodes = linspace(0,L,T.n_nodes)';
T.elements = [(1:(T.n_nodes-1))', (2:T.n_nodes)' ]; 
T.edges = (1:T.n_edges)'; 

T.edgeFlags = zeros(T.n_edges,1);
T.edgeFlags(1) = 1;
T.edgeFlags(end) =1;

