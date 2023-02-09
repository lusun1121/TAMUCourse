function [local_rhs] = local_assemblyrhs(localnodes,f,t,FE_at_Quad, Quad, p)

mat_B = [(localnodes(2,:) - localnodes(1,:))',(localnodes(3,:) - localnodes(1,:))']; %[(dim)x(dim)]
det_B = det(mat_B);  

q_points = (mat_B*Quad.xhat' + repmat(localnodes(1,:)',1,Quad.nq))'; 
rhs_vals=f(q_points(:,1),q_points(:,2),t*ones(Quad.nq,1));

local_rhs = zeros(3,1);

for q_index = 1:Quad.nq % run through quadrature points on element

    local_rhs = local_rhs + rhs_vals(q_index)...              
                          * FE_at_Quad.hat_phi(q_index,:)' ... 
                          * Quad.what(q_index)...              
                          * abs(det_B);                            

end