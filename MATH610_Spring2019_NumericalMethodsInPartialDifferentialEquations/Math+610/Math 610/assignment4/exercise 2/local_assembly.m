function [cell_matrix] = local_assembly(localnodes,FE_at_Quad, Quad)

mat_B = [(localnodes(2,:) - localnodes(1,:))',(localnodes(3,:) - localnodes(1,:))']; %[(dim)x(dim)]
det_B = det(mat_B);  
inv_B = 1/det_B*[mat_B(2,2),-mat_B(1,2);-mat_B(2,1),mat_B(1,1)];

cell_matrix = zeros(3,3);

for q_index = 1:Quad.nq 
    
    grad_phi_at_q_point = [FE_at_Quad.hat_phix(q_index,:)',...
                            FE_at_Quad.hat_phiy(q_index,:)']*inv_B;
                        
    grad_phi_ij_matrix = grad_phi_at_q_point*grad_phi_at_q_point';
    
    cell_matrix = cell_matrix + grad_phi_ij_matrix...
                              * Quad.what(q_index) ...
                              * abs(det_B);
                          
end

