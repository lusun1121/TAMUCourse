%  Written by Spencer Patty
%  Jan 20, 2017
%  for use in Math 610 Lab
function [local_rhs] = assembleLocalRhs(f, localnodes, FE_at_Quad, Quad)
% Given an element defined by nodes in localnodes, we assemble the local rhs
% term associated with this element.
%
% In 1D, the element is an interval, so there are p=2 nodes on the interval
% and p=2 linear basis functions corresponding to those nodes
%
% Input:
% f: @(x) rhs function
% localnodes: [pxdim] is a column vector of vertices on element (each vertex is a "row" vector)
% Quad.nq: [1x1] number of quadrature points on reference element
% Quad.what: [nqx1] weights on reference element for each quad point
% Quad.xhat: [nqxdim] coordinates of quad point on reference element
% FE_at_Quad.hat_phi: [nqxp] p basis functions on reference element
%                            evaluated at the Quad.nq quad points
% FE_at_Quad.hat_phix: [nqxp] x derivative of p basis functions on 
%                             reference element evaluated at the Quad.nq 
%                             quad points
%
% Output:
% local_rhs = [px1] rhs contributions from all interactions of basis function
%
% Let localnodes = [p1; p2 ;...; pp] where each pi [1xdim] is the ith node
% on local cell T
%
% In 1D, localnodes = [p1; p2] = [2x1]
% In 2D, localnodes = [p1; p2; p3]  = [3x2], etc

% 1D Conversion from reference element, That = (0,1) to T = (p1, p2)
% (we use points are row vectors)
%
% x = F(xhat) = xhat*B' + c = (p2-p1)*xhat + p1
%
% Note that 
%
% uhat(xhat) := u(F(xhat)) so that grad_uhat = B' grad_u  or grad_u = B'\grad_uhat
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
    
% precompute various things at the quadrature points if needed
q_points = Quad.xhat*mat_B' + repmat(localnodes(1,:),Quad.nq,1); % [nqx2] list of quadrature points in cell 
rhs_vals = f(q_points);

% preallocate space for cell_matrix
local_rhs = zeros(2,1);

for q_index = 1:Quad.nq % run through quadrature points on element
    
    % exploiting matlab's speed with vector operations, we avoid doing
    % the for loops to construct the local rhs.  Below in
    % the commented section is the same code using for loops.

    local_rhs = local_rhs + rhs_vals(q_index)...              % f(q_point)
                        * FE_at_Quad.hat_phi(q_index,:)'... % bases on reference element at q_point [2x1]
                        * Quad.what(q_index)...             % quadrature weight
                        * det_B;                            % reference area-element transform   


%     % The following shows the local rhs with the for loop.
%     for j = 1:2
%         local_rhs(j) = local_rhs(j) + ...
%             (rhs_vals(q_index)*FE_at_Quad.hat_phi(q_index,j))*Quad.what(q_index)*det_B;
%     end
%     

end
