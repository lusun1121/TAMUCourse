function [rhs_g] = assembleLocalRhsg1(g, v1,v2, FE_at_Quad, Quad)

rhs_g = zeros(2,1);

for q_index = 1:Quad.nq 

    rhs_g = rhs_g + g([Quad.xhat(q_index);Quad.xhat(q_index)].*(v2'-v1')+v1')...               
                          * FE_at_Quad.hat_phi(q_index,:)' ... 
                          * Quad.what(q_index)...             
                          * norm(v2-v1);                             
end