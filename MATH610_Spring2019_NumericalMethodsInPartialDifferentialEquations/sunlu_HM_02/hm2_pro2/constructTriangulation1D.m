%  Written by Spencer Patty
%  Jan 20, 2017
%  for use in Math 610 Lab
function [T] = constructTriangulation1D(L, num_elements,num_start,num_end)
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

% Setup edges flags:  (recall in 1D an edge is a vertex node)
%   0 = non boundary edge
%   1 = Dirichlet edge
%   2 = Natural edge
%   3 = Robin edge
T.edgeFlags = zeros(T.n_edges,1);
T.edgeFlags(1) = num_start; % set Dirichlet
T.edgeFlags(end) = num_end; % set Dirichlet
T.meshsize = L/T.n_elements;
%T.edgeFlags(1) = 1; % set Dirichlet
%T.edgeFlags(end) = 1; % set Dirichlet
T.dirichletdofs = find(T.edgeFlags == 1); % using the edgeFlags  
                                                   % scheme and assuming 1D,
                                                   % so edges and nodes have
                                                   % the same ordering for 
                                                   % vertices and no interior 
                                                   % nodes can be
                                                   % Dirichlet.
%DoFHandler.naturaldofs = find(T.edgeFlags == 2); % using the edgeFlags                                                     
% Extract the coordinates of the dirichlet dofs so that we can easily
% apply the dirichlet bc to them.  In 1D they are only in vertices, but in 
% higher dimensions, we would have to construct all the edge dofs
% pertaining to the dirichlet edges.
T.dirichletdofs_coordinates = T.nodes(T.dirichletdofs,:);
%DoFHandler.naturaldofs_coordinates = T.nodes(DoFHandler.naturaldofs,:);                                                   
% enumerate the Free DoFs = all dofs that are not Dirichlet

T.freedofs = [find(T.edgeFlags ==0); ((T.n_nodes+1):T.n_nodes)'];


