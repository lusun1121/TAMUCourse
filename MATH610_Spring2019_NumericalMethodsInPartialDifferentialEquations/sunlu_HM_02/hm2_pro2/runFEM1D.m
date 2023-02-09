function [method, mesh_size, n_dofs, Linferror,uh] = runFEM1D(num_elements,method,epsillion) 
u_exact = @(x)(1-exp(-x/epsillion))/(1-exp(-1/epsillion));
uh = FD_pro2(epsillion, num_elements,method);
L = 1;


mesh_size = L/num_elements;
n_dofs = num_elements + 1;

n_nodes = mesh_size .* (0:n_dofs-1);


local_Linf_error = 0;

for q_index = 1: n_dofs % run through quadrature points on element
    diff_val = uh(q_index) - u_exact(n_nodes(q_index)); %[1x1] 
    local_Linf_error = max(local_Linf_error, abs(diff_val));
end
Linferror = local_Linf_error;

switch method
    case 1
        method_str ='unpwind';
    case 2
        method_str ='downwind';
    case 3
        method_str ='centered';
end
if (epsillion ==1)
    ep='epmax';
else
    ep='epmin';
end
    


%i=1:1:num_elements-1;
%nodes= (i-1)*mesh_size;
%lmbd = figure;
%plot(nodes,uh(i),'+',nodes,u_exact(nodes),'-');
%xlabel('x');
%ylabel('u(x)');
%legend('finite difference approximations','exact solution');
%title({[method_str],['epsilion = ', num2str(epsillion)],['number of intervals =',num2str(num_elements)]});
%saveas(lmbd, ['\\blender\homes\l\u\lusun8825\nt\AccountSettings\Desktop\610\FEM1D_For_Students\sunlu\result\',...
%    'EX2_hm2_',ep,'_mesh',num2str(num_elements),'_',method_str,'.png'],'png');
end