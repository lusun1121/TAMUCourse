u_exact = @(x)sin(pi*x(1))*sin(pi*x(2)); 
grad_u_exact = @(x) [pi*cos(pi*x(1))*sin(pi*x(2)),...
                     pi*sin(pi*x(1))*cos(pi*x(2))]; 
f = @(x)4*pi^2*sin(pi*x(1))*sin(pi*x(2));  
g_D = @(x) u_exact(x);
g=@(x)(x(1)==-1)*2*pi*sin(pi*x(2))+(x(1)==0)*2*pi*sin(pi*x(2))-(x(1)==1)*2*pi*sin(pi*x(2))...
    -(x(2)==1)*2*pi*sin(pi*x(1))-(x(2)==0)*2*pi*sin(pi*x(1))+(x(2)==-1)*2*pi*sin(pi*x(1));
T=gettriangle(1);

p = 1; 
DoFHandler = constructDoFHandler(T,p);

uh = zeros(DoFHandler.n_dofs,1);
RHS = zeros(DoFHandler.n_dofs,1);
A = spalloc(DoFHandler.n_dofs, DoFHandler.n_dofs, 15*T.n_nodes);

quad_n_points = 4;
Quad = getQuadOnRefElement(quad_n_points);

[ FE_at_Quad] = feEval( Quad, p );


for cell = 1:T.n_elements
    dofIndices = DoFHandler.dofs(cell,:); 
    vertices = T.nodes(T.elements(cell,:),:);
    
    [cell_matrix,cell_rhs,~] = local_assembly(vertices,f, FE_at_Quad, Quad, p); 
    
    A(dofIndices,dofIndices) = A(dofIndices,dofIndices) + cell_matrix;
    RHS(dofIndices) = RHS(dofIndices) + cell_rhs;
end

Quad1 = getQuadOnRefElement1D(quad_n_points);
[FE_at_Quad1] = feEval1D( Quad1, p );

for cell=1:T.n_boundaryedges
    dofIndices=T.edges(T.boundaryedges(cell),:);
    v1=T.nodes(T.edges(T.boundaryedges(cell),1),:);
    v2=T.nodes(T.edges(T.boundaryedges(cell),2),:);
    [rhs_g] = assembleLocalRhsg1(g, v1,v2, FE_at_Quad1, Quad1);
     RHS(dofIndices) = RHS(dofIndices) +rhs_g;
end


uh(2:end) = A(2:end, 2:end) \ RHS (2:end);