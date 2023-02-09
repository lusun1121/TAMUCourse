%  Written by Spencer Patty
%  Jan 20, 2017
%  for use in Math 610 Lab
function [local_Linf_error] = computeLocalInfErrors(localnodes, uh_at_localnodes, u_exact, Quad, FE_at_Quad)
% Given an element defined by nodes in localnodes, we compute the Linfinity
% error compared to the exact solutions on the local element.
%
% In 1D, the element is an interval, so there are p=2 nodes on the interval
% and 2 basis functions
%
% Input:
% localnodes: [p x dim] is a column vector of nodes on element (each node is a "row" vector)
% uh_at_localnodes: [px1] vector of values of uh at p localnodes
% u_exact: @(x) exact function
% Quad : a reference element Quadrature for calculating the L inf norm
% FE_at_Quad : the Finite element ref shape functions evaluated at Quad
%
%
% Output:

% local_Linf_error: [1x1]
%
% Let localnodes = [p1; p2 ;...; pp] where each pi [1xdim] is the ith node
% on local cell T
%
% In 1D, localnodes = [p1; p2] = [2x1]
% In 2D, localnodes = [p1; p2; p3]  = [3x2], etc

% 1D Conversion from reference element, That = (0,1) to T = (p1, p2)
% (we use points are row vectors)
%
% x = F(xhat) = xhat*B + c = (p2-p1)*xhat + p1
%
% 1D Basis functions on reference element:
% (we use gradients are column vectors)
%
% lambdahat_1 = x          grad lambdahat_1 = 1
% lambdahat_2 = 1-x        grad lamddahat_2 = -1
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
mat_B = (localnodes(2,:) - localnodes(1,:))'; %[(p-1)x(p-1)]
det_B = abs(det(mat_B));    
q_points = Quad.xhat*mat_B' + repmat(localnodes(1,:),Quad.nq,1); % [nqx(p-1)] list of quadrature points in cell 

% evaluate solutions and gradients at quadrature points
u_exact_at_q_points = u_exact(q_points); %[nqx1]
uh_at_q_points = FE_at_Quad.hat_phi*uh_at_localnodes;           %[nqx1] = [nqxp] [px1]

% assemble error analysis
local_Linf_error = 0;

for q_index = 1:Quad.nq % run through quadrature points on element
    diff_val = uh_at_q_points(q_index) - u_exact_at_q_points(q_index); %[1x1] 
    local_Linf_error = max(local_Linf_error, abs(diff_val));
end




