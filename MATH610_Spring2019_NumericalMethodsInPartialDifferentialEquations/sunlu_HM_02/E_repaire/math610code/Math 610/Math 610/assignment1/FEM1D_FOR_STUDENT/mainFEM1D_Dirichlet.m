
L = 2;

u_exact = @(x)sin(pi*(x(1)+x(2)));       % if x input is [nxdim] returns [nx1]
grad_u_exact = @(x) [pi*cos(pi*(x(1)+x(2)));...
                     pi*cos(pi*(x(1)+x(2)))]; % if x input is [nxdim] returns [nxdim]

% right hand side function
f = @(x)u_exact(x);  % = - u''(x)   if x input is [nxdim]  returns [nx1]

% begin FEM code

T = TriangleReader('LShapedDomain.3.node','LShapedDomain.3.ele','LShapedDomain.3.edge');

%
% setup_system
%
p = 1;  % polynomial degree of lagrange finite element basis
DoFHandler = constructDoFHandler(T,p);

uh = zeros(DoFHandler.n_dofs,1);
RHS = zeros(DoFHandler.n_dofs,1);
A = spalloc(DoFHandler.n_dofs, DoFHandler.n_dofs, 15*T.n_nodes); % upper bound on NNZ in matrix for 1D is 2*p+1 interactions per node (p on either side plus itself)

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
    
    [cell_matrix,cell_rhs,~] = local_assamble(vertices, FE_at_Quad, Quad, p); 
    
    % contribute local terms to global stiffness and RHS structures
    A(dofIndices,dofIndices) = A(dofIndices,dofIndices) + cell_matrix;
    RHS(dofIndices) = RHS(dofIndices) + cell_rhs;
end

% solve_system
%
% Only solve system for unconstrained dofs (the non Dirichlet ones)
uh = A \ RHS;

%compute error
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




