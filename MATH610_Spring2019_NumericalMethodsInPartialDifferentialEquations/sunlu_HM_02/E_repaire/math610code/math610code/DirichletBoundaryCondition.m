
L2=zeros(4,1);
H1=zeros(4,1);
Linf=zeros(4,1);
E=zeros(4,6);

u_exact = @(x,y)cos(2*pi*x).*cos(2*pi*y); 
grad_u_exact = @(x,y) [-2*pi*sin(2*pi*x).*cos(2*pi*y),...
                     -2*pi*cos(2*pi*x).*sin(2*pi*y)]; 

%u_exact_0 = @(r,theta)r^(2/3).*sin(2/3*theta);
%u_exact = @(x,y)u_exact((x^2+y^2)^(1/2),atan2(y,x));
%grad_u_exact_0 = @(r,theta)[2/3*r^(-1/3).*sin(2/3*theta),2/3*r^(2/3).*cos(2/3*theta)]*[cos(theta),sin(theta);-1/r,1/r];
%grad_u_exact=@(x,y)grad_u_exact_0((x^2+y^2)^(1/2),atan2(y,x));

% right hand side function
f = @(x,y)(8*pi^2)*(cos(2*pi*x).*cos(2*pi*y));  
g_D = @(x,y) u_exact(x,y);
%f_0 = @(r,theta) -2/9*r^(-4/3).*sin(2/3*theta) + 1/r.*2/3*r^(-1/3).*sin(2/3*theta)+1/r^2.*(-4/9)*r^(2/3).*sin(2/3*theta);%u/x''+u/y''=u/r''+1/r*u/r'+1/r^2*u/theta''
%f = @(x,y)f_0((x^2+y^2)^(1/2),atan2(y,x));


for j=1:4
    T=gettriangle(j);

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
quad_n_points = 4;%1;%2;%3;%5??
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

uh(DoFHandler.dirichletdofs) = g_D(DoFHandler.dirichletdofs_coordinates(:,1),DoFHandler.dirichletdofs_coordinates(:,2));

RHS = RHS - A*uh;


%
% solve_system
%
% Only solve system for unconstrained dofs (the non Dirichlet ones)
uh(DoFHandler.freedofs) = A(DoFHandler.freedofs, DoFHandler.freedofs) \ RHS (DoFHandler.freedofs);

%compute error
Quad_Error = getQuadOnRefElement(quad_n_points);
FE_at_Quad_Error = feEval(Quad_Error, p);

% In the L inf error, we don't need a weight, only a bunch of xhat locations
% so we combine the Quad_Error xhats with a uniformly distributed set of
% nodes through the reference element to make a nice L inf quadrature rule


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
    
    localLinf = computeLocalInfErrors(vertices, uh(dofIndices), u_exact, Quad_Error, FE_at_Quad_Error,p);
    Linferror = max( Linferror, localLinf); 
end
L2(j) = sqrt(L2sqrd);
H1(j) = sqrt(H1sqrd);
Linf(j)=Linferror;
fprintf('%1.4e %1.4e %1.4e\n',double(Linf(j)),double(L2(j)),double(H1(j)));
end

h=[0.2,0.1,0.05,0.025];
for i=1:4
    E(i,:)=[Linf(i),0,L2(i),0,H1(i),0];
    
   if (i>1)
	E(i,2) = (E(i,1)-E(i-1,1))/(h(i)-h(i-1)); %Linf error
       E(i,4) = (E(i,3)-E(i-1,3))/(h(i)-h(i-1)); %L2 error
       E(i,6) = (E(i,5)-E(i-1,5))/(h(i)-h(i-1)); %H1 error

       %E(i,2) = log(E(i,1)/E(i-1,1))/log(h(i)/h(i-1)); %Linf error
       %E(i,4) = log(E(i,3)/E(i-1,3))/log(h(i)/h(i-1)); %L2 error
       %E(i,6) = log(E(i,5)/E(i-1,5))/log(h(i)/h(i-1)); %H1 error
   end
    

   fprintf('%f %f %f %f %f %f\n', E(i,1),E(i,2),E(i,3),E(i,4),E(i,5),E(i,6));


if j ==4
	plotfunction;
end


end

loglog(h,E(:,1),h,E(:,3),h,E(:,5));
legend('Inf','L2','H1');

