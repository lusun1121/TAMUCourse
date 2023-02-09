function [mass_matrix,stiffness_matrix] = local_assembly(localnodes,FE_at_Quad, Quad, p)

mat_B = [(localnodes(2,:) - localnodes(1,:))',(localnodes(3,:) - localnodes(1,:))']; %[(dim)x(dim)]
det_B = det(mat_B);  
inv_B = mat_B\eye(2);

% preallocate space for cell_matrix
mass_matrix = zeros(3,3);
stiffness_matrix = zeros(3,3);

for q_index = 1:Quad.nq % run through quadrature points on element
    
    grad_phi_at_q_point = [FE_at_Quad.hat_phix(q_index,:)',...
                            FE_at_Quad.hat_phiy(q_index,:)']*inv_B;%3x2
                        
    grad_phi_ij_matrix = grad_phi_at_q_point*grad_phi_at_q_point';%3x3
    
    stiffness_matrix = stiffness_matrix + grad_phi_ij_matrix...
                              * Quad.what(q_index) ...
                              * abs(det_B);
                          
    phi_ij_matrix =  FE_at_Quad.hat_phi(q_index,:)'*...
                     FE_at_Quad.hat_phi(q_index,:); %3x3                 
                          
    mass_matrix = mass_matrix + phi_ij_matrix...
                            * Quad.what(q_index) ...
                            * abs(det_B);

end

