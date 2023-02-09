function [mesh_size, n_dofs, L2error, H1error, Linferror] = runFEM1Dproblem(n_elements,p,quad_degree1,quad_degree2,bPlotSolutions)

L = 50;

q=200;
S=100;
D=8.8*10^7;
a=S*L^2/D;
Q=q*L^2/(2*D);



u_exact=@(x)Q/a*(1-1/sinh(sqrt(a))*(sinh(sqrt(a)*x/L)+sinh(sqrt(a)*(1-x/L))));


grad_u_exact = @(x)Q/a*(-1/sinh(sqrt(a))*(sqrt(a)/L*(cosh(sqrt(a)*x/L)-...
    cosh(sqrt(a)*(1-x/L)))));

f=@(x)q/(2*D)*ones(size(x));

g_D = @(x) u_exact(x); 



T = constructTriangulation1D(L, n_elements);


DoFHandler = constructDoFHandler(T,p);

uh = zeros(DoFHandler.n_dofs,1);
RHS = zeros(DoFHandler.n_dofs,1);
A1 = spalloc(DoFHandler.n_dofs, DoFHandler.n_dofs, (2*p+1)*T.n_nodes);
A2 = spalloc(DoFHandler.n_dofs, DoFHandler.n_dofs, (2*p+1)*T.n_nodes);

m=max(quad_degree1,quad_degree2);
Quad1 = getQuadOnRefElement(quad_degree1);
Quad2 = getQuadOnRefElement(quad_degree2);
Quadm = getQuadOnRefElement(m);


[ FE_at_Quad1] = feEval( Quad1, p );
[ FE_at_Quad2] = feEval( Quad2, p );
[ FE_at_Quadm] = feEval( Quadm, p );



for cell = 1:T.n_elements
    dofIndices = DoFHandler.dofs(cell,:);
    vertices = T.nodes(T.elements(cell,:),:); 
    
    cell_matrix1 = assembleLocalStiffness(vertices, FE_at_Quad2, Quad2, p);
    cell_matrix2 = assembleLocalmass(vertices, FE_at_Quad1, Quad1, p);
    cell_rhs    = assembleLocalRhs(f, vertices, FE_at_Quadm, Quadm,p); 
    
   
    A1(dofIndices,dofIndices) = A1(dofIndices,dofIndices) + cell_matrix1;
    A2(dofIndices,dofIndices) = A2(dofIndices,dofIndices) + cell_matrix2;
    RHS(dofIndices) = RHS(dofIndices) + cell_rhs;
end

A=A1+S/D*A2;

uh(DoFHandler.dirichletdofs) = g_D(DoFHandler.dirichletdofs_coordinates);
RHS = RHS - A*uh;



uh(DoFHandler.freedofs) = A(DoFHandler.freedofs, DoFHandler.freedofs) \ RHS (DoFHandler.freedofs);

Quad_Error = getQuadOnRefElement(m);
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
    
    localLinf = computeLocalInfErrors(vertices, uh(dofIndices), u_exact, QuadInf_Error, FE_at_QuadInf_Error);
    Linferror = max( Linferror, localLinf); 
end
L2error = sqrt(L2sqrd);
H1error = sqrt(H1sqrd);

mesh_size = (L-0)/T.n_elements;
n_dofs = DoFHandler.n_dofs;

% fprintf('The Linf error for N = %d dofs, h=%1.4e, is %1.4e.\n',DoFHandler.n_dofs,mesh_size, Linferror);
% fprintf('The L2   error for N = %d dofs, h=%1.4e, is %1.4e.\n',DoFHandler.n_dofs,mesh_size, L2error);
% fprintf('The H1   error for N = %d dofs, h=%1.4e, is %1.4e.\n',DoFHandler.n_dofs,mesh_size, H1error);

%
% output_solution
%
if (bPlotSolutions)
    figure;
    plot( T.nodes, uh(1:T.n_nodes), T.nodes, u_exact(T.nodes) ); % only plot the vertices with linear solution
    xlabel('x');
    ylabel('solution');
    title('1D Laplacian with Dirichlet BC');
    legend('uh', 'uexact');

%     figure(2);
%     % plot errors at half nodes since likely we are very close a true nodes, so
%     % we plot at half nodes to get a better glimpse of what the error looks
%     % like.
%     nodes_half = 0.5*(T.nodes(2:T.n_nodes) + T.nodes(1:T.n_nodes-1) );
%     uh_half = 0.5*( uh(2:T.n_nodes) + uh(1:T.n_nodes-1) );
% 
%     plot( nodes_half, uh_half - u_exact(nodes_half) );
%     xlabel('x');
%     ylabel('uh-u');
%     title('Error at half nodes');
end