%  Written by Spencer Patty
%  Jan 20, 2017
%  for use in Math 610 Lab
function [local_rhs] = assembleLocalRhs(f, localnodes, FE_at_Quad, Quad, p)
% Given an element defined by nodes in localnodes, we assemble the local rhs
% term associated with this element.
%
% In 1D, the element is an interval, so there are (dim+1)= 2 nodes on the interval
% and p+1 Lagrange basis functions of degree p on each interval
%
% Input:
% f: @(x) rhs function
% localnodes: [(dim+1)xdim] is a column vector of vertices on element (each vertex is a "row" vector)
% Quad.nq: [1x1] number of quadrature points on reference element
% Quad.what: [nqx1] weights on reference element for each quad point
% Quad.xhat: [nqxdim] coordinates of quad point on reference element
% FE_at_Quad.hat_phi: [nqx(p+1)] p+1 basis functions on reference element
%                                evaluated at the Quad.nq quad points
% FE_at_Quad.hat_phix: [nqx(p+1)] x derivative of p+1 basis functions on 
%                                 reference element evaluated at the Quad.nq 
%                                 quad points
%
% Output:
% local_rhs = [(p+1)x1] rhs contributions from all interactions of basis function
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
% Thus the grad dot grad term from integration by parts of laplacian 
% transformed onto integration over the reference element is
%
% \int_{T} f(x) lambda_i dx
%   = \int_{That} f(F(xhat)) lambdahat_i(xhat)  |det(B)| dxhat
%


%affine map info
mat_B = (localnodes(2,:) - localnodes(1,:))';
det_B = abs(det(mat_B));
% inv_B = 1/mat_B;

% precompute various things at the quadrature points if needed
q_points = Quad.xhat*mat_B + repmat(localnodes(1,:),Quad.nq,1); % [nqxdim] list of quadrature points in cell 
rhs_vals = f(q_points);

% preallocate space for cell_matrix
local_rhs = zeros(p+1,1);

for q_index = 1:Quad.nq % run through quadrature points on element
    
    % exploiting matlab's speed with vector operations, we avoid doing
    % the for loops to construct the local rhs.  Below in
    % the commented section is the same code using for loops.

    local_rhs = local_rhs + rhs_vals(q_index)...               % f(q_point)
                          * FE_at_Quad.hat_phi(q_index,:)' ... % bases on reference element at q_point [(p+1)x1]
                          * Quad.what(q_index)...              % quadrature weight
                          * det_B;                             % reference area-element transform   


%     % The following shows the local rhs with the for loop.
%     for j = 1:(p+1)
%         local_rhs(j) = local_rhs(j) + ...
%             (rhs_vals(q_index)*FE_at_Quad.hat_phi(q_index,j))*Quad.what(q_index)*det_B;
%     end
%     

end
