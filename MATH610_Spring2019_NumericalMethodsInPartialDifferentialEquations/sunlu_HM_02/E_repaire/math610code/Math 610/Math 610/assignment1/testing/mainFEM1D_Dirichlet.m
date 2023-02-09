%  Written by Spencer Patty
%  Jan 20, 2017
%  for use in Math 610 Lab
%
%  Constructs a 1D linear FEM to solve 
%  the Laplacian problem
%
%   -u''(x) = f(x)   for x in (0,L)
%   
%   with Dirichlet BC:
%
%    u(0) = a  = g_D(0)
%    u(L) = b  = g_D(L)
%
%  Weak Formulation of PDE:
%  Let \mathcal{T}_h := [0=x_0, x_1 , ... , x_N=L]  be the triangulation.
%  Lambda := (0,L)
%  partialLambda := {0,L} = boundary of domain
%  partialLambdaD := {x \in partialLambda with dirichlet bc being applied }
%  Let Vh  := { vh \in \mathbb{P}^1(\mathcal{T}_h) | vh(xD) = g_D(xD) for all xD  in partialLambdaD}
%      Vh0 := { vh \in \mathbb{P}^1(\mathcal{T}_h) | vh(xD) = 0 for all xD in partialLambdaD }
%
%  Find uh in Vh such that for all vh in Vh0
%
%     a(uh,vh) = L(vh)
%
%  where
%
%  a(uh,vh) = sum_{T in \mathcal{T}_h}  \int_{T} grad uh . grad vh dx     
%
%  L(vh) = sum_{T in \mathcal{T}_h} \int_{T} f(x) vh(x) dx
%
%
% Language of elements:
%  1D element is an interval and 1D boundary element is a node
%  2D element is a triangle and 2D boundary element is an edge
%  3D element is a tetrahedron and 3D boundary element is a face
%

% Define problem domain
L = 2;

% exact solution to compare our solution with in processing
u_exact = @(x) sin(4*x)/2;       % if x input is [nxdim] returns [nx1]
grad_u_exact = @(x) 2*cos(4*x);  % if x input is [nxdim] returns [nxdim]

% right hand side function
f = @(x) 8*sin(4*x);  % = - u''(x)   if x input is [nxdim]  returns [nx1]

% Dirichlet boundary conditions function
g_D = @(x) u_exact(x);  % if x input is [nxdim] returns [nx1]


%
% begin FEM code
%

%
% triangulate domain in 1D
%
num_elements = 800; %100;
T = constructTriangulation1D(L, num_elements);

%
% setup_system
%
p = 1;  % polynomial degree of lagrange finite element basis
DoFHandler = constructDoFHandler(T,p);

uh = zeros(DoFHandler.n_dofs,1);
RHS = zeros(DoFHandler.n_dofs,1);
A = spalloc(DoFHandler.n_dofs, DoFHandler.n_dofs, (2*p+1)*T.n_nodes); % upper bound on NNZ in matrix for 1D is 2*p+1 interactions per node (p on either side plus itself)

% Choose quadrature on reference element and evaluate finite element shape 
% functions on reference element at quadrature points
quad_n_points = 4;
Quad = getQuadOnRefElement(quad_n_points);

% Evaluate shape value and shape gradient(dim components) on reference 
% element at quadrature points for the bulk elements.
% each term in FE_at_Quad structure is [nqx1]
[ FE_at_Quad] = feEval( Quad, p );

% If needed for stiffness or RHS, define a quadrature rule for boundary 
% element or edges of elements and call it QuadEdge that can be passed to 
% assembly functions.  Then evaluate the FE_at_QuadEdge ... for a reference
% boundary element or edge element.


%
% assemble_system
%

% Loop through all cells and construct the local stiffness and rhs
% quantities, then distribute them to the global matrix and rhs vectors.
for cell = 1:T.n_elements
    dofIndices = DoFHandler.dofs(cell,:); % [1x(p+1)]  extract indices pertaining to cell nodes
    vertices = T.nodes(T.elements(cell,:),:); % [(dim+1)xdim] coordinates of vertices
    
    cell_matrix = assembleLocalStiffness(vertices, FE_at_Quad, Quad, p); % assemble local stiffness terms
    cell_rhs    = assembleLocalRhs(f, vertices, FE_at_Quad, Quad,p); % assemble local rhs terms
    
    % contribute local terms to global stiffness and RHS structures
    A(dofIndices,dofIndices) = A(dofIndices,dofIndices) + cell_matrix;
    RHS(dofIndices) = RHS(dofIndices) + cell_rhs;
end


%
% apply_boundary_conditions
%

% Apply Dirichlet BC by setting solution values
% then updating the RHS vector to reflect that information
% so we can solve the restricted problem on freeNodes properly
uh(DoFHandler.dirichletdofs) = g_D(DoFHandler.dirichletdofs_coordinates);
RHS = RHS - A*uh;


%
% solve_system
%
% Only solve system for unconstrained dofs (the non Dirichlet ones)
uh(DoFHandler.freedofs) = A(DoFHandler.freedofs, DoFHandler.freedofs) \ RHS (DoFHandler.freedofs);
% uh(DoFHandler.freedofs) = pcg(A(DoFHandler.freedofs, DoFHandler.freedofs) , RHS (DoFHandler.freedofs), 1e-14, 4*DoFHandler.n_dofs);

%
% process_system and compute errors
%

% Choose a possibly different quadrature for error analysis and evaluate
% the fe shape functions at this new quadrature.  Typically we want one
% that is different and exact for higher dof so we don't just sample at the
% places where the solution is most accurate, but get a good representation
% of the solution through the domain.
Quad_Error = getQuadOnRefElement(quad_n_points);
FE_at_Quad_Error = feEval(Quad_Error, p);

% In the L inf error, we don't need a weight, only a bunch of xhat locations
% so we combine the Quad_Error xhats with a uniformly distributed set of
% nodes through the reference element to make a nice L inf quadrature rule
n_inf_nodes = 10;
QuadInf_Error.nq = Quad_Error.nq+n_inf_nodes;
QuadInf_Error.xhat = [Quad_Error.xhat; linspace(0,1,n_inf_nodes)'];
FE_at_QuadInf_Error = feEval(QuadInf_Error, p);


% Loop through the cells and compute the local errors and aggregate them
% according to the norms
L2sqrd = 0;
H1sqrd = 0;
Linferror = 0;
for cell = 1:T.n_elements
    dofIndices = DoFHandler.dofs(cell,:); % [1x(p+1)]  extract indices pertaining to cell nodes
    vertices = T.nodes(T.elements(cell,:),:); % [(dim+1)xdim] coordinates of vertices
    
    [localL2sqrd, localH1sqrd] = computeLocalErrors(vertices, uh(dofIndices), u_exact, grad_u_exact, Quad_Error, FE_at_Quad_Error, p);
    L2sqrd = L2sqrd + localL2sqrd;
    H1sqrd = H1sqrd + localH1sqrd;
    
    localLinf = computeLocalInfErrors(vertices, uh(dofIndices), u_exact, QuadInf_Error, FE_at_QuadInf_Error,p);
    Linferror = max( Linferror, localLinf); 
end
L2error = sqrt(L2sqrd);
H1error = sqrt(H1sqrd);

mesh_size = (L-0)/T.n_elements;

fprintf('The Linf error for N = %d dofs, h=%1.4e, is %1.4e.\n',DoFHandler.n_dofs,mesh_size,double(Linferror));
fprintf('The L2   error for N = %d dofs, h=%1.4e, is %1.4e.\n',DoFHandler.n_dofs,mesh_size, double(L2error));
fprintf('The H1   error for N = %d dofs, h=%1.4e, is %1.4e.\n',DoFHandler.n_dofs,mesh_size, double(H1error));

%
% output_solution
%
figure(1);
plot( T.nodes, uh(1:T.n_nodes), T.nodes, u_exact(T.nodes) ); % only plot the vertices with linear solution
xlabel('x');
ylabel('solution');
title('1D Laplacian with Dirichlet BC');
legend('uh', 'uexact');


figure(2);
% plot errors at half nodes since likely we are very close a true nodes, so
% we plot at half nodes to get a better glimpse of what the error looks
% like.
nodes_half = 0.5*(T.nodes(2:T.n_nodes) + T.nodes(1:T.n_nodes-1) );
uh_half = 0.5*( uh(2:T.n_nodes) + uh(1:T.n_nodes-1) );

plot( nodes_half, uh_half - u_exact(nodes_half) );
xlabel('x');
ylabel('uh-u');
title('Error at half nodes');




