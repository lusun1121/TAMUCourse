function [rhs_g] = assembleLocalRhsg(g, localnodes, FE_at_Quad, Quad)

mat_B = (localnodes(2,:) - localnodes(1,:))';
det_B = abs(det(mat_B));

q_points = Quad.xhat*mat_B + repmat(localnodes(1,:),Quad.nq,1); 
rhs_vals = g(q_points);

rhs_g = zeros(2,1);

for q_index = 1:Quad.nq 

    rhs_g = rhs_g + rhs_vals(q_index)...               
                          * FE_at_Quad.hat_phi(q_index,:)' ... 
                          * Quad.what(q_index)...              
                          * det_B;                             
end