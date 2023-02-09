
L = pi;


u_exact = @(x) cos(2*x)+1;       
grad_u_exact = @(x) -2*sin(2*x);  


f = @(x) 7*cos(2*x)+3;  


g_D = @(x) u_exact(x); 


num_elements = 800;
T = constructTriangulation1D(L, num_elements);


p = 1;  
DoFHandler = constructDoFHandler(T,p);

uh = zeros(DoFHandler.n_dofs,1);
RHS = zeros(DoFHandler.n_dofs,1);
A1 = spalloc(DoFHandler.n_dofs, DoFHandler.n_dofs, (2*p+1)*T.n_nodes); 
A2 = spalloc(DoFHandler.n_dofs, DoFHandler.n_dofs, (2*p+1)*T.n_nodes); 


quad_n_points = 2;
Quad = getQuadOnRefElement(quad_n_points);


[ FE_at_Quad] = feEval( Quad, p );

for cell = 1:T.n_elements
    dofIndices = DoFHandler.dofs(cell,:); 
    vertices = T.nodes(T.elements(cell,:),:); 
    
    cell_matrix1 = assembleLocalStiffness(vertices, FE_at_Quad, Quad, p); 
    cell_rhs    = assembleLocalRhs(f, vertices, FE_at_Quad, Quad,p); 
    cell_matrix2 = assembleLocalmass(vertices, FE_at_Quad, Quad, p);
    
    A1(dofIndices,dofIndices) = A1(dofIndices,dofIndices) + cell_matrix1;
    A2(dofIndices,dofIndices) = A2(dofIndices,dofIndices) + cell_matrix2;
    RHS(dofIndices) = RHS(dofIndices) + cell_rhs;
end
A=A1+3*A2;


uh(DoFHandler.dirichletdofs) = g_D(DoFHandler.dirichletdofs_coordinates);
RHS = RHS - A*uh;



uh(DoFHandler.freedofs) = A(DoFHandler.freedofs, DoFHandler.freedofs) \ RHS (DoFHandler.freedofs);

Quad_Error = getQuadOnRefElement(quad_n_points);
FE_at_Quad_Error = feEval(Quad_Error, p);


n_inf_nodes = 10;
QuadInf_Error.nq = Quad_Error.nq+n_inf_nodes;
QuadInf_Error.xhat = [Quad_Error.xhat; linspace(0,1,n_inf_nodes)'];
FE_at_QuadInf_Error = feEval(QuadInf_Error, p);



L2sqrd = 0;
H1sqrd = 0;
Linferror = 0;
for cell = 1:T.n_elements
    dofIndices = DoFHandler.dofs(cell,:); 
    vertices = T.nodes(T.elements(cell,:),:); 
    
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

figure(1);
plot( T.nodes, uh(1:T.n_nodes), T.nodes, u_exact(T.nodes) ); % only plot the vertices with linear solution
xlabel('x');
ylabel('solution');
title('1D Laplacian with Dirichlet BC');
legend('uh', 'uexact');


figure(2);

nodes_half = 0.5*(T.nodes(2:T.n_nodes) + T.nodes(1:T.n_nodes-1) );
uh_half = 0.5*( uh(2:T.n_nodes) + uh(1:T.n_nodes-1) );

plot( nodes_half, uh_half - u_exact(nodes_half) );
xlabel('x');
ylabel('uh-u');
title('Error at half nodes');
