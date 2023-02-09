
f = @(x,y,t)ones(size(x));%nx1 
%get triangle and degree of freedom
for j=1:3
    
T=gettriangle(j);
p = 1;  % polynomial degree of lagrange finite element basis
DoFHandler = constructDoFHandler(T,p);

M = spalloc(DoFHandler.n_dofs, DoFHandler.n_dofs, 15*T.n_nodes);
S = spalloc(DoFHandler.n_dofs, DoFHandler.n_dofs, 15*T.n_nodes);

quad_n_points = 4;
Quad = getQuadOnRefElement(quad_n_points);
[FE_at_Quad] = feEval( Quad, p );

% assemble_system
for cell = 1:T.n_elements
    dofIndices = DoFHandler.dofs(cell,:); % [1x(p+1)]  extract indices pertaining to cell nodes
    vertices = T.nodes(T.elements(cell,:),:); % [(dim+1)xdim] coordinates of vertices
    
    [mass_matrix,stiffness_matrix] = local_assembly(vertices, FE_at_Quad, Quad, p); 
    
    % contribute local terms to global stiffness and RHS structures
    M(dofIndices,dofIndices) = M(dofIndices,dofIndices) + mass_matrix;
    S(dofIndices,dofIndices) = S(dofIndices,dofIndices) + stiffness_matrix;
end

t=linspace(0,3,79);
tau=3/78;
uh1 = zeros(DoFHandler.n_dofs,1);
uh2 =  zeros(DoFHandler.n_dofs,1);


%backward Euler
for i=1:79
    ti=t(i);
    RHS = zeros(DoFHandler.n_dofs,1);
    for cell = 1:T.n_elements
    dofIndices = DoFHandler.dofs(cell,:); 
    vertices = T.nodes(T.elements(cell,:),:); 
    local_rhs=local_assemblyrhs(vertices,f,ti,FE_at_Quad, Quad, p);
    RHS(dofIndices)=RHS(dofIndices)+local_rhs;
    end
    
    uh2=(1-tau)*uh1-tau*(M\S)*uh1+tau*M\RHS;
    
    %compute error
    if (ti==1 || ti==3)
        fprintf('plot the function\n');
        figure
        plotfunction;
    end
    
    uh1=uh2;
end

end