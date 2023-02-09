
%Input the element
L2error=zeros(4,1);
for i=1:4
T = gettriangle(i);

%u=@(x)abs(x(1)-x(2));
%u=@(x)sin(pi*(x(1)+x(2)));
u=@(x)double(x(1)^2+x(2)^2>=1/4);
quad_n_points=3;
L2sqrd = 0;
p=1;

for ele=1:T.n_elements
    %vertical [v1;v2;v3]
    v = T.nodes(T.elements(ele,:),:);
    %The coefficients for basis function
    uh=[u(v(1,:));u(v(2,:));u(v(3,:))];
    Quad_Error = getQuadOnRefElement(quad_n_points);
    FE_at_Quad_Error = feEval(Quad_Error, p);
    [localL2sqrd] = compute_error(v, uh, u, Quad_Error, FE_at_Quad_Error, p);
    L2sqrd = L2sqrd + localL2sqrd;
end
L2error(i) = sqrt(L2sqrd);
fprintf('%f\n',L2error(i));
end
h=[0.2,0.1,0.05,0.025];
for i=2:4
    rate=(L2error(i)-L2error(i-1))/(h(i)-h(i-1));%rate=log(L2error(i)/L2error(i-1))/log(h(i)/h(i-1));
    fprintf('%f\n',rate);
end