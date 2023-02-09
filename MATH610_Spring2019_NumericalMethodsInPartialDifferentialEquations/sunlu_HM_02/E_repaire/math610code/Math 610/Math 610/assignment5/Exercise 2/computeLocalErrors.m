function [local_L2_error_sqrd, local_H1_error_sqrd] = computeLocalErrors(localvertices, uh_at_localdofs, u_exact, grad_u_exact,t, Quad, FE_at_Quad, p)

mat_B = [(localvertices(2,:) - localvertices(1,:))',(localvertices(3,:) - localvertices(1,:))']; %[(dim)x(dim)]
det_B = det(mat_B); 
inv_B = 1/det_B*[mat_B(2,2),-mat_B(1,2);-mat_B(2,1),mat_B(1,1)];

q_points = (mat_B*Quad.xhat' + repmat(localvertices(1,:)',1,Quad.nq))'; % [nqxdim] 

u_exact_at_q_points=u_exact(q_points(:,1),q_points(:,2),t*ones(Quad.nq,1));
grad_u_exact_at_q_points=grad_u_exact(q_points(:,1),q_points(:,2),t*ones(Quad.nq,1));

uh_at_q_points = FE_at_Quad.hat_phi*uh_at_localdofs; %[nqx1] 

gradx_uh_at_qhat_points =[FE_at_Quad.hat_phix*uh_at_localdofs,...
                          FE_at_Quad.hat_phiy*uh_at_localdofs];
                      
grad_uh_at_q_points = gradx_uh_at_qhat_points*inv_B; %[nqx2] 

% assemble error analysis
local_L2_error_sqrd = 0;
local_H1semi_error_sqrd = 0;

for q_index = 1:Quad.nq

    diff_val = uh_at_q_points(q_index) - u_exact_at_q_points(q_index); %[1x1]
    diff_grad = grad_uh_at_q_points(q_index,:) - grad_u_exact_at_q_points(q_index,:); % [1xdim]
    
    local_L2_error_sqrd = local_L2_error_sqrd + diff_val*diff_val * Quad.what(q_index)*abs(det_B);
    local_H1semi_error_sqrd = local_H1semi_error_sqrd + dot(diff_grad,diff_grad)*Quad.what(q_index)*abs(det_B);
    
end
local_H1_error_sqrd = local_H1semi_error_sqrd + local_L2_error_sqrd;


