
function [cell_matrix, area] = assembleLocalStiffness2(localnodes, FE_at_Quad, Quad)

mat_B = (localnodes(2,:) - localnodes(1,:))';
det_B = abs(det(mat_B));
    

cell_matrix = zeros(2,2);
area = det_B;

for q_index = 1:Quad.nq 
    
    grad_phi_at_q_point=[FE_at_Quad.hat_phi(q_index,:)];
  

    grad_phi_ij_matrix = grad_phi_at_q_point'*grad_phi_at_q_point;  

    cell_matrix = cell_matrix + grad_phi_ij_matrix' * Quad.what(q_index) * det_B;   
  
end
