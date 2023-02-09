function [cell_matrix, area] = assembleLocalmasstiffness(localnodes, FE_at_Quad, Quad,p)

mat_B = (localnodes(2,:) - localnodes(1,:))';
det_B = abs(det(mat_B));
inv_B = 1/mat_B;

cell_matrix = zeros(p+1,p+1);
area = det_B;

for q_index = 1:Quad.nq 
    
    phi_at_q_point=[FE_at_Quad.hat_phi(q_index,:)];
    grad_phi_at_q_point = inv_B*FE_at_Quad.hat_phix(q_index,:);
  

    phi_ij_matrix =  grad_phi_at_q_point'*phi_at_q_point;  

    cell_matrix = cell_matrix + phi_ij_matrix' * Quad.what(q_index) * det_B;   
  
end