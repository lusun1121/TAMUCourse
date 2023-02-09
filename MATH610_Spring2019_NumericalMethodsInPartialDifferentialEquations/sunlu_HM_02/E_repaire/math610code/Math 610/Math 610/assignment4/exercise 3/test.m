
f = @(x)1;
g_D = @(x)0;
estimate=zeros(4,1);
for j=1:4
%T=gettriangle1(j);
T=gettriangle2(j);
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
    
    [cell_matrix,cell_rhs,~] = local_assembly(vertices,f, FE_at_Quad, Quad, p); 
    
    % contribute local terms to global stiffness and RHS structures
    A(dofIndices,dofIndices) = A(dofIndices,dofIndices) + cell_matrix;
    RHS(dofIndices) = RHS(dofIndices) + cell_rhs;
end

% solve_system
% 
% Only solve system for unconstrained dofs (the non Dirichlet ones)

for t=1:T.n_boundarynodes
uh(DoFHandler.dirichletdofs(t)) = 0;
end

RHS = RHS - A*uh;


%
% solve_system
%
% Only solve system for unconstrained dofs (the non Dirichlet ones)
uh(DoFHandler.freedofs) = A(DoFHandler.freedofs, DoFHandler.freedofs) \ RHS (DoFHandler.freedofs);

%compute eta
epsilon=0;
for ele=1:T.n_elements
    vertices = T.nodes(T.elements(ele,:),:);
    h_T=dot(vertices(1,:)-vertices(2,:),vertices(1,:)-vertices(2,:))^(1/2);
    mat_B = [(vertices(2,:) - vertices(1,:))',...
             (vertices(3,:) - vertices(1,:))']; 
    %f_norm =h_T*sqrt(abs(det(mat_B))/2);
    f_norm=abs(det(mat_B));
    jump=0;
for e=1:3
    edge=abs(T.elementToEdgeMap(ele,e));
    v=T.nodes(T.edges(edge,:)',:);
    h=dot(v(1,:)-v(2,:),v(1,:)-v(2,:))^(1/2);
    dofIndices=abs(T.edgesToElementsMap(edge,:));
    if dofIndices(2)==0
        grad_uh_T1=0;
        grad_uh_T2=0;
    else
    grad_uh_T1=grad_at_element(T.nodes(T.elements(dofIndices(1),:),:),...
                                          uh(T.elements(dofIndices(1),:)));
    grad_uh_T2=grad_at_element(T.nodes(T.elements(dofIndices(2),:),:),...
                                          uh(T.elements(dofIndices(2),:)));
    end
    jump=jump+h^(1/2)*(dot(grad_uh_T1-grad_uh_T2,grad_uh_T1-grad_uh_T2))^(1/2);
end
    eta=f_norm+jump;
    epsilon=epsilon+eta^2;
end
estimate(j)=sqrt(epsilon);
end

%compute rate
rate=zeros(4,1);
for i=2:4
    rate(i)=log(estimate(i)/estimate(i-1))/log(1/2);
end
