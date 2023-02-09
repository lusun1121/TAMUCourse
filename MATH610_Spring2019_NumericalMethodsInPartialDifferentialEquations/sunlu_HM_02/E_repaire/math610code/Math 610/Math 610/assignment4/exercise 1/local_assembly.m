function [cell_matrix, local_rhs,area] = local_assembly(localnodes,f,FE_at_Quad, Quad, p)

mat_B = [(localnodes(2,:) - localnodes(1,:))',(localnodes(3,:) - localnodes(1,:))']; %[(dim)x(dim)]
det_B = det(mat_B);  
inv_B = 1/det_B*[mat_B(2,2),-mat_B(1,2);-mat_B(2,1),mat_B(1,1)];

q_points = (mat_B*Quad.xhat' + repmat(localnodes(1,:)',1,Quad.nq))'; % [nqxdim] list of real quadrature points in cell not ref element
rhs_vals=zeros(Quad.nq,1);
for t=1:Quad.nq
rhs_vals(t) = f(q_points(t,:));
end

% preallocate space for cell_matrix
cell_matrix = zeros(3,3);
local_rhs = zeros(3,1);
area = abs(det_B);

for q_index = 1:Quad.nq % run through quadrature points on element
    
    grad_phi_at_q_point = [FE_at_Quad.hat_phix1(q_index,:)*inv_B;...
                           FE_at_Quad.hat_phix2(q_index,:)*inv_B;...
                           FE_at_Quad.hat_phix3(q_index,:)*inv_B];
    
    grad_phi_ij_matrix = grad_phi_at_q_point*grad_phi_at_q_point';  
   
    cell_matrix = cell_matrix +2*grad_phi_ij_matrix...
                              * Quad.what(q_index) ...
                              * abs(det_B);

    local_rhs = local_rhs + rhs_vals(q_index)...               
                          * FE_at_Quad.hat_phi(q_index,:)' ... 
                          * Quad.what(q_index)...              
                          * abs(det_B);                           

end

