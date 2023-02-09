function [local_L2_error_sqrd] = compute_error(localvertices, uh_at_localdofs, u_exact, Quad, FE_at_Quad, p)
% Given an element defined by nodes in localnodes, we compute the L2 and H1
% errors compared to the exact solutions on the local element.

% Input:
% localvertices: [(dim+1) x dim] is a column vector of vertices on element (each node is a "row" vector)
% uh_at_localdofs: [(p+1)*(p+2)/2)] vector of values of uh at p+dim degrees of
%                            freedom
% u_exact: @(x) exact function
% grad_u_exact: @(x) gradient of exact function
% Quad : a reference element Quadrature for calculating the L2, H1 norms
% FE_at_Quad : the Finite element ref shape functions evaluated at Quad
%
% Output:
% local_L2_error_sqrd: [1x1] 
% local_H1_error_sqrd: [1x1]
%
% vertex on local cell T
%
% In 1D, localnodes = [v_1; v_2] = [2x1]
% In 2D, localnodes = [v_1; v_2; v_3]  = [3x2], etc

mat_B = [(localvertices(2,:) - localvertices(1,:));(localvertices(3,:) - localvertices(1,:))]; %[(dim)x(dim)]
det_B = det(mat_B);    

q_points = Quad.xhat*mat_B + repmat(localvertices(1,:),Quad.nq,1); % [nqxdim] list of real quadrature points in cell not ref element

% evaluate solutions and gradients at quadrature points
u_exact_at_q_points=zeros(Quad.nq,1);

for t=1:Quad.nq
u_exact_at_q_points(t) = u_exact(q_points(t,:)); %[nqx1]
end

uh_at_q_points = FE_at_Quad.hat_phi*uh_at_localdofs; %[nqx1]
% assemble error analysis
local_L2_error_sqrd = 0;

for q_index = 1:Quad.nq

    diff_val = uh_at_q_points(q_index) - u_exact_at_q_points(q_index); %[1x1]
    
    local_L2_error_sqrd = local_L2_error_sqrd + diff_val*diff_val * Quad.what(q_index)*det_B;
    
end