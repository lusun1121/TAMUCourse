function [mesh_size, n_dofs, L2error, H1error, Linferror] = runFEM1D(num_elements,p,quad_n_points, bPlotSolutions,problem_case,S)

%  Constructs a 1D linear FEM to solve 
%  the Laplacian problem
%
%   K(x)*u''(x)+Q(x)*u(x) = f(x)   for x in (0,L)
%   
%   with Dirichlet BC:
%
%    u(0) = a  = g_D(0)
%    u(L) = b  = g_D(L)
%
%   with mixed BC( Dirichlet + natural):
%
%    u(0) = a = g_D(0)
%    u'(L)= 0
%
%
%  Weak Formulation of PDE:
%  Let \mathcal{T}_h := [0=x_0, x_1 , ... , x_N=L]  be the triangulation.
%  Lambda := (0,L)
%  partialLambda := {0,L} = boundary of domain
%  partialLambdaD := {x \in partialLambda with dirichlet bc being applied }
%  Let Vh  := { vh \in \mathbb{P}^1(\mathcal{T}_h) | 
%                           vh(xD) = g_D(xD) for all xD  in partialLambdaD}
%      Vh0 := { vh \in \mathbb{P}^1(\mathcal{T}_h) | 
%                            vh(xD) = 0 for all xD in partialLambdaD }
%
%  Find uh in Vh such that for all vh in Vh0
%
%     a(uh,vh) = L(vh)
%
%  where
%
%  a(uh,vh) =  sum_{T in \mathcal{T}_h} \int_{T} (-K) * grad uh . grad vh dx 
%            + sum_{T in \mathcal{T}_h} \int_{T} (+Q) *  uh .  vh dx
%
%  L(vh) = sum_{T in \mathcal{T}_h} \int_{T} f(x) vh(x) dx
%
%
% Language of elements:
%  1D element is an interval and 1D boundary element is a node
%  2D element is a triangle and 2D boundary element is an edge
%  3D element is a tetrahedron and 3D boundary element is a face
%
%
%num_elements : 25 or 50 or 100 or 200 (number of intervals)
%choose       : [problem_case, S]
%
%problem_case   problem
% 0             homework4
% 1             homework1
% 2             homework2
% 3             homework3.1
% 4             homework3.2 
%
%S : 100 or 1000 or 10000 (parameter)
%
%
%problem_case=choose(1);
%S=choose(2);

% Define problem domain and parameter

if(problem_case~=0)
    L = 50;
    q = 200;
    D = 8.8*10^7;
    a = S*L^2/D;
    b = q*L^4/2/D;
    Q = q*L^2/2/D;
    
    function_K = @(x)-1;        %K(x)
    function_Q = @(x)S/D;       %Q(x)

    
    switch(problem_case)
        
        case 1
            % exact solution to compare our solution with in processing
            u_exact = @(x)b/a*(-(x/L).^2+(x/L)-2/a+2/(a*sinh(sqrt(a)))...
                .*(sinh(sqrt(a)*(x/L))+sinh(sqrt(a)*(1-x/L))));
            % if x input is [nxdim] returns [nx1]
            %  right hand side function
            f = @(x)q/2/D*x.*(L-x);% if x input is [nxdim]  returns [nx1]
            num_start = 1;     %'1': Dirichelet
            num_end = 1; 
            
        case 2
            u_exact = @(x)Q/a*(1-1/sinh(sqrt(a))*(sinh(sqrt(a)*x/L)+...
                sinh(sqrt(a)*(1-x/L))));
            f = @(x)q/2/D*x./x;
            num_start = 1;     
            num_end = 1;    
            
        case 3
            u_exact = @(x)b/a*(-(x/L).^2+(x/L)-2/a+1/(a*cosh(sqrt(a)))...
                .*(sqrt(a)*sinh(sqrt(a)*(x/L))+2*cosh(sqrt(a)*(1-x/L))));
            f = @(x)q/2/D*x.*(L-x);
            num_start = 1;     
            num_end = 0;       
            %'2': Natural but has same function as '0': Unbounded
            
        case 4
            u_exact = @(x)Q/a*(1+sinh(sqrt(a))/cosh(sqrt(a))*...
                sinh(sqrt(a)*x/L)-cosh(sqrt(a)*x/L));
            f = @(x)q/2/D*x./x;
            num_start = 1;     
            num_end = 0;  
            
    end
    
    syms t;
    grad_u_exact = @(x)vpa(subs(diff(u_exact(t)),t,x));
    % if x input is [nxdim] returns [nxdim]
    
    % Dirichlet boundary conditions function
    g_D = @(x)0;% if x input is [nxdim] returns [nx1]
    
else
    L=1;
    function_K   = @(x)(x>=0).*(x<pi/6).*1 + (x>pi/6).*(x<pi/4).*2 + ...
        (x>pi/4).*(x<=1).*3;
    function_Q   = @(x)0;
    u_exact      = @(x)(x>=0).*(x<pi/6).*(12/pi).*x +...
        (x>pi/6).*(x<pi/4).*((6/pi)*x+1) + (x>pi/4).*(x<=1).*((4/pi)*x+3/2);
    grad_u_exact = @(x)(x>=0).*(x<pi/6).*(12/pi) +...
        (x>pi/6).*(x<pi/4).*(6/pi) + (x>pi/4).*(x<=1).*(4/pi);
    
    f   = @(x)0*x./x;
    g_D =@(x)u_exact(x);
    
    num_start = 1;     %Dirichelet
    num_end   = 1;     %Dirichelet
end
%
% begin FEM code
%
%
% triangulate domain in 1D
%
T = constructTriangulation1D(L, num_elements,num_start,num_end);
%
% setup_system
%
%p = 1;  % polynomial degree of lagrange finite element basis
DoFHandler = constructDoFHandler(T,p);

uh = zeros(DoFHandler.n_dofs,1);
RHS = zeros(DoFHandler.n_dofs,1);
A = spalloc(DoFHandler.n_dofs, DoFHandler.n_dofs, (2*p+1)*T.n_nodes); 
% upper bound on NNZ in matrix for 1D is 2*p+1 interactions per node 
%(p on either side plus itself)

% Choose quadrature on reference element and evaluate finite element shape 
% functions on reference element at quadrature points
%quad_n_points = 4;
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
    %When K(x) is not constant, function_x = ( x_i + x_(i+1) ) / 2
    %i.e.mid point of the interval I_ei = [ x_i, x_(i+1) ]
    %K_ei: = K(function_x); namely, mid_point rule
    function_x=L/T.n_elements*(cell-1/2); 
    
    dofIndices = DoFHandler.dofs(cell,:); 
    % [1x(p+1)]  extract indices pertaining to cell nodes
    
    vertices = T.nodes(T.elements(cell,:),:); 
    % [(dim+1)xdim] coordinates of vertices
    
    cell_matrix1 = assembleLocalStiffness(vertices, FE_at_Quad, Quad, p);
    % assemble local stiffness terms
    
    cell_matrix = -function_K(function_x)*cell_matrix1.fi1 ...
                  +function_Q(function_x)*cell_matrix1.fi0;
    
    cell_rhs    = assembleLocalRhs(f, vertices, FE_at_Quad, Quad,p); 
    % assemble local rhs terms
    
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
uh(DoFHandler.freedofs) = A(DoFHandler.freedofs, DoFHandler.freedofs) ...
                          \ RHS (DoFHandler.freedofs);
% uh(DoFHandler.freedofs) = pcg(A(DoFHandler.freedofs, DoFHandler.freedofs) 
%                  ,RHS (DoFHandler.freedofs), 1e-14, 4*DoFHandler.n_dofs);

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
    dofIndices = DoFHandler.dofs(cell,:); 
    % [1x(p+1)]  extract indices pertaining to cell nodes
    vertices = T.nodes(T.elements(cell,:),:); 
    % [(dim+1)xdim] coordinates of vertices
    
    if (problem_case==0)
        localL2sqrd = computeLocalErrors(vertices, uh(dofIndices), ...
            u_exact, grad_u_exact, Quad_Error, FE_at_Quad_Error, p);
    else
        [localL2sqrd, localH1sqrd] = computeLocalErrorsH1(vertices, ...
            uh(dofIndices), u_exact, grad_u_exact, Quad_Error, FE_at_Quad_Error, p);
        H1sqrd = H1sqrd + localH1sqrd;
    end
    L2sqrd = L2sqrd + localL2sqrd;

    
    localLinf = computeLocalInfErrors(vertices, uh(dofIndices), ...
        u_exact, QuadInf_Error, FE_at_QuadInf_Error,p);
    Linferror = max( Linferror, localLinf); 
end
L2error = sqrt(L2sqrd);
%
% output_solution
%
if(problem_case~=0)
    H1error = sqrt(H1sqrd);
else
    H1error = false;
end

mesh_size = (L-0)/T.n_elements;
n_dofs = DoFHandler.n_dofs;


%
% output_solution
%
if (bPlotSolutions)
    figure(1);
    plot( T.nodes, uh(1:T.n_nodes), T.nodes, u_exact(T.nodes) ); 
    % only plot the vertices with linear solution
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
end




