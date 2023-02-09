L2=zeros(4,1);
H1=zeros(4,1);
Linf=zeros(4,1);
E=zeros(4,6);

u_exact = @(x)sin(pi*x(1))*sin(pi*x(2)); 
grad_u_exact = @(x) [pi*cos(pi*x(1))*sin(pi*x(2)),...
                     pi*sin(pi*x(1))*cos(pi*x(2))]; 
f = @(x)4*pi^2*sin(pi*x(1))*sin(pi*x(2));  
g_D = @(x) u_exact(x);

for j=1:4

T=gettriangle(j);

p = 1;  % polynomial degree of lagrange finite element basis
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
%     [a,b,k,l,g]=getfandinterval(v1,v2,dofIndices(1),dofIndices(2));
    [a,b,g]=getfandintervaltext(v1,v2);
    localnodes=[a;b];
    [rhs_g] = assembleLocalRhsg(g, localnodes, FE_at_Quad1, Quad1);
     RHS(dofIndices) = RHS(dofIndices) +rhs_g;
end


uh(2:end) = A(2:end, 2:end) \ RHS (2:end);


%compute the mean value and get the constant
mean=0;
for cell = 1:T.n_elements
    dofIndices = DoFHandler.dofs(cell,:); 
    vertices = T.nodes(T.elements(cell,:),:);
    mean=mean+meanatelement(vertices,uh(dofIndices),Quad,FE_at_Quad);
end

uh=uh+(4/pi^2-mean)/3;

%compute error
Quad_Error = getQuadOnRefElement(quad_n_points);
FE_at_Quad_Error = feEval(Quad_Error, p);
QuadInf_Error.nq = Quad_Error.nq+10;
addpoints=[0,0;1/3,0;2/3,0;1,0;0,1/3;1/3,1/3;2/3,1/3;0,2/3;1/3,2/3;0,1];
QuadInf_Error.xhat = [Quad_Error.xhat;addpoints];
FE_at_QuadInf_Error = feEval(QuadInf_Error, p);

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
    
    localLinf = computeLocalInfErrors(vertices, uh(dofIndices), u_exact, QuadInf_Error, FE_at_QuadInf_Error,p);
    Linferror = max( Linferror, localLinf); 
end
L2(j) = sqrt(L2sqrd);
H1(j) = sqrt(H1sqrd);
Linf(j)=Linferror;
uh(1)
fprintf('%1.4e %1.4e %1.4e\n',double(Linf(j)),double(L2(j)),double(H1(j)));

end

h=[0.2,0.1,0.05,0.025];
for i=1:4
    E(i,:)=[Linf(i),0,L2(i),0,H1(i),0];
    
   if (i>1)
       E(i,2) = log(E(i,1)/E(i-1,1))/log(h(i)/h(i-1)); %Linf error
       E(i,4) = log(E(i,3)/E(i-1,3))/log(h(i)/h(i-1)); %L2 error
       E(i,6) = log(E(i,5)/E(i-1,5))/log(h(i)/h(i-1)); %H1 error
   end
    

   fprintf('%f %f %f %f %f %f\n', E(i,1),E(i,2),E(i,3),E(i,4),E(i,5),E(i,6));
end

loglog(h,E(:,1),h,E(:,3),h,E(:,5));
legend('Inf','L2','H1');
