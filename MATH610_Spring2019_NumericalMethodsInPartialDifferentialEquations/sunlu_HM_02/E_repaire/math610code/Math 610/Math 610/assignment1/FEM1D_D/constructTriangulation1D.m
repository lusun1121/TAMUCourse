%  Written by Spencer Patty
%  Jan 20, 2017
%  for use in Math 610 Lab
function [T] = constructTriangulation1D(L, num_elements)
%
% We create a structure, T which houses all the parts of the triangulation
% as say T.n_elements  or T.n_edges  or T.nodes... That way the
% triangulation is completely inclosed and can easily be passed around if
% needed
%
%  The domain is (0,L) and internally here, we choose where the Dirichlet
%  and Natural BC will be applied
%

T.n_elements = num_elements;
T.n_nodes = T.n_elements+1;
T.n_edges = T.n_nodes;
T.nodes = linspace(0,L,T.n_nodes)'; % [n_nodes x dim]
T.elements = [(1:(T.n_nodes-1))', (2:T.n_nodes)' ]; % [n_elements x dim+1 ]
T.edges = (1:T.n_edges)'; % [n_edges x dim]

% Setup boundary element flags:
%   0 = non boundary edge
%   1 = Dirichlet edge
%   2 = Natural edge
%   3 = Robin edge
T.bdElementFlags = zeros(T.n_elements+1,1);
T.bdElementFlags(1) = 1; % set Dirichlet
T.bdElementFlags(end) = 1; % set Dirichlet

% pointers to different classes of Nodes
T.boundaryNodes = [1; T.n_nodes]; % [n_boundaryNodes x 1]
T.interiorNodes = (2:(T.n_nodes-1))'; % [n_interiorNodes x 1]
T.boundaryElementDirichlet = T.boundaryNodes; % [ n_boundaryElementDirichlet x 1] both boundary Elements are Dirichlet
T.boundaryElementNatural = [];  % [n_boundaryElementNatural x 1] no Boundary Elements have Natural bc applied

% extract the unique Nodes on the Dirichlet Boundary Elements
% In 1D, the boundary elements are nodes so there is nothing to do.  In 
% higher dimensions, we extract the nodes and then reduce to a unique
% set of node pointers to nodes array
T.DirichletNodes = T.boundaryElementDirichlet;

% Extract the set of non-Dirichlet nodes on which we will solve the system
% after boundary conditions have been set = interiorNodesPtr + non dirichlet
% boundary nodes
T.freeNodes = T.interiorNodes;  