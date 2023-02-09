function [DoFHandler] = constructDoFHandler(T, p)
%
% Input:
% T = Triangulation structure 
% p = polynomial degree on each element
%
% Output:
% DoFHandler.p = polynomial degree of finite element basis on ref element
% DoFHandler.dofs = [T.n_elements x p+1 ] for each element, the columns are
%                   the indices of dofs in the global system
% DoFHandler.n_dofs = number of global dofs
% DoFHandler.dirichletdofs = list of indices of dofs where Dirichlet bc
%                            will be applied
% DoFHandler.dirichletdof_coordinates = coordinates of dofs where dirichlet
%                                       bc will be applied
% DoFHandler.freedofs = list of indices of unconstrained dofs 
%
% 

%
% The global ordering of dofs is done by multiple passes through the
% elements enumerating the dofs from each set as we go:
%
% In 1D, we would order them 1) vertices, then 2) interior of elements
%
% In 2D, we would order them 1) vertices, then 2) dofs on edges, then 
%   3) dofs on interiors of elements
%
% In 3D, we would order them 1) vertices, then 2) dofs on faces (on edges 
% of faces then interior of face), then 3) dofs on interior of element
%
%
%

% store polynomial degree of lagrange fe basis on each element
DoFHandler.p = p;

%
% enumerate the dofs
%
DoFHandler.dofs = zeros(T.n_elements, (p+1)*(p+2)/2);

% first the existing vertices
DoFHandler.dofs(:,1:3) = T.elements; 

% add the interior dofs for each cell
count = T.n_nodes+1;
% add loop through each element here to add the extra p-1 dofs on interior
% of each cell.


% total number of degrees of freedom (constrained and unconstrained
DoFHandler.n_dofs = count-1;


%
% Construct lists of constrained and un constrained dofs
%

% enumerate the Dirichlet DoFs

% T.edgeFlags scheme:
%   0 = non boundary edge
%   1 = Dirichlet edge
%   2 = Natural edge
%   3 = Robin edge

% DoFHandler.dirichletdofs = [1; T.n_nodes]; %  first and last vertex (dirichlet on both ends)
DoFHandler.dirichletdofs = T.boundarynodes(1);
                                                   
% Extract the coordinates of the dirichlet dofs so that we can easily
% apply the dirichlet bc to them.  In 1D they are only in vertices, but in 
% higher dimensions, we would have to construct all the edge dofs
% pertaining to the dirichlet edges.
DoFHandler.dirichletdofs_coordinates = T.nodes(DoFHandler.dirichletdofs,:);
                                                   
%enumerate the Free DoFs = all dofs that are not Dirichlet
if (p>1)
    DoFHandler.freedofs = [find(T.edgeFlags ~=1); ((T.n_nodes+1):DoFHandler.n_dofs)'];
else
    DoFHandler.freedofs = [2:T.n_nodes]';
end

