u_exact = @(x,y,t)t.*exp(-t).*cos(3*pi*x).*cos(pi*y);%nx1
grad_u_exact = @(x,y,t) [-3*pi*t.*exp(-t).*sin(3*pi*x).*cos(pi*y),...
                         -pi*t.*exp(-t).*cos(3*pi*x).*sin(pi*y)];%nx2
% right hand side function
f = @(x,y,t)(1+4*t+10*pi^2*t).*exp(-t).*cos(3*pi*x).*cos(pi*y);%nx1 
g_D = @(x,y,t) u_exact(x,y,t);
%get triangle and degree of freedom
output=zeros(3,8);
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

t=linspace(0,3,19);
tau=3/18;
uh1 = zeros(DoFHandler.n_dofs,1);
uh2 =  zeros(DoFHandler.n_dofs,1);


%backward Euler
for i=2:19
    ti=t(i);
    RHS = zeros(DoFHandler.n_dofs,1);
    for cell = 1:T.n_elements
    dofIndices = DoFHandler.dofs(cell,:); 
    vertices = T.nodes(T.elements(cell,:),:); 
    local_rhs=local_assemblyrhs(vertices,f,ti,FE_at_Quad, Quad, p);
    RHS(dofIndices)=RHS(dofIndices)+local_rhs;
    end
    
    uh2=((1+5*tau)*M+tau*S)\(M*uh1+tau*RHS);
    
    %compute error
    if (ti==1 || ti==3)
        Quad_Error = getQuadOnRefElement(quad_n_points);
        FE_at_Quad_Error = feEval(Quad_Error, p);
        
        L2sqrd = 0;
        H1sqrd = 0;
     for cell = 1:T.n_elements
     dofIndices = DoFHandler.dofs(cell,:); 
     vertices = T.nodes(T.elements(cell,:),:); 
    
     [localL2sqrd, localH1sqrd] = computeLocalErrors(vertices, uh2(dofIndices), u_exact, grad_u_exact,ti,Quad_Error, FE_at_Quad_Error, p);
     L2sqrd = L2sqrd + localL2sqrd;
     H1sqrd = H1sqrd + localH1sqrd;
     end
     L2 = sqrt(L2sqrd);
     H1 = sqrt(H1sqrd);
     if ti==1
     output(j,1)=L2;
     output(j,5)=H1;
     fprintf('lalalalalala\n');
     else
     output(j,3)=L2;
     output(j,7)=H1;
     end
    end
    
    uh1=uh2;
end

end

for i=2:3
    for j=2:2:8
        output(i,j)=log(output(i,j-1)/output(i-1,j-1))/log(1/2);
    end
end





