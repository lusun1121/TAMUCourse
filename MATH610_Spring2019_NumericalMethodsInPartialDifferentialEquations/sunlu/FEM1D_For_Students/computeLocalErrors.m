%  Written by Spencer Patty
%  Jan 20, 2017
%  for use in Math 610 Lab
function [local_L2_error_sqrd] = computeLocalErrors(localvertices, uh_at_localdofs, u_exact, grad_u_exact, Quad, FE_at_Quad, p)
%function [local_L2_error_sqrd, local_H1_error_sqrd] = computeLocalErrors(localvertices, uh_at_localdofs, u_exact, grad_u_exact, Quad, FE_at_Quad, p)
% Given an element defined by nodes in localnodes, we compute the L2 and H1
% errors compared to the exact solutions on the local element.
%
% In 1D, the element is an interval, so there are (dim+1)=2 nodes on the interval
% and p+1 lagrange basis functions of degree p
%
% Input:
% localvertices: [(dim+1) x dim] is a column vector of vertices on element (each node is a "row" vector)
% uh_at_localdofs: [(p+1)x1] vector of values of uh at p+1 degrees of
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
% Let localnodes = [v_1; ...; v_{dim+1}] where each v_i [1xdim] is the ith
% vertex on local cell T
%
% In 1D, localnodes = [v_1; v_2] = [2x1]
% In 2D, localnodes = [v_1; v_2; v_3]  = [3x2], etc

% 1D Conversion from reference element, That = (0,1) to T = (v_1, v_2)
% (we use points are row vectors)
%
% x = F(xhat) = xhat*B + c = (v_2-v_1)*xhat + v_1
%
% Note that 
%
% uhat(xhat) := u(F(xhat)) so that grad_uhat = B grad_u  or grad_u = B\grad_uhat
%  
% and
%
% dx = |det(B)| dxhat
%


%affine map info
mat_B = (localvertices(2,:) - localvertices(1,:))'; %[(dim)x(dim)]
det_B = abs(det(mat_B));    
%inv_B = 1/mat_B;

q_points = Quad.xhat*mat_B + repmat(localvertices(1,:),Quad.nq,1); % [nqxdim] list of real quadrature points in cell not ref element

% evaluate solutions and gradients at quadrature points
u_exact_at_q_points = u_exact(q_points); %[nqx1]
%grad_u_exact_at_q_points = grad_u_exact(q_points); %[nqxdim]

uh_at_q_points = FE_at_Quad.hat_phi*uh_at_localdofs; %[nqx1] = [nqx(p+1)][(p+1)x1]

% Assemble components of gradient on reference element
%gradx_uh_at_qhat_points = FE_at_Quad.hat_phix*uh_at_localdofs; %[nqx1] = [nqx(p+1)][(p+1)x1]

% combine gradient components into gradient and transform onto real element
%grad_uh_at_q_points = (inv_B*([gradx_uh_at_qhat_points]'))'; %[nqxdim] = ([dimxdim] [dimxnq] )'

% assemble error analysis
local_L2_error_sqrd = 0;
%local_H1semi_error_sqrd = 0;

for q_index = 1:Quad.nq

    diff_val = uh_at_q_points(q_index) - u_exact_at_q_points(q_index); %[1x1]
    %diff_grad = grad_uh_at_q_points(q_index) - grad_u_exact_at_q_points(q_index); % [1xdim]
    
    local_L2_error_sqrd = local_L2_error_sqrd + diff_val*diff_val * Quad.what(q_index)*det_B;
    %local_H1semi_error_sqrd = local_H1semi_error_sqrd + dot(diff_grad,diff_grad)*Quad.what(q_index)*det_B;
    
end
%local_H1_error_sqrd = local_H1semi_error_sqrd + local_L2_error_sqrd;



