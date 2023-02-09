function uh = FD_pro2(epsillion, num_elements,method)
%epsillion = 1; %=10^(-2);
b = -1/epsillion;
L = 1;
num_start = 1;
num_end = 1;
%num_elements = 10;% 20; 40; 80; 160; 320; 640;
T = constructTriangulation1D( L, num_elements, num_start, num_end);

u_exact = @(x)(1-exp(-x/epsillion))/(1-exp(-1/epsillion));

f = @(x)0.*x;
g_D = @(x) x;

T = constructTriangulation1D( L, num_elements, num_start, num_end);

uh = zeros(T.n_nodes,1);
RHS = zeros(T.n_nodes,1);
A = spalloc(T.n_nodes,T.n_nodes, 3*T.n_nodes); 

A_matrix = Differnece_One(A,T,method,b);



for cell = 1:(T.n_nodes)
    RHS = f(T.nodes(cell));
end

uh(T.dirichletdofs) = g_D(T.dirichletdofs_coordinates);
RHS = RHS - A_matrix*uh;
uh(T.freedofs) = A_matrix(T.freedofs, T.freedofs) ...
                          \ RHS (T.freedofs);
                      



