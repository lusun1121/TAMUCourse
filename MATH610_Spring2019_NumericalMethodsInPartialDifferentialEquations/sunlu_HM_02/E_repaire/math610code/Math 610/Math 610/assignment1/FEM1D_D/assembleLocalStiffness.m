%  Written by Spencer Patty
%  Jan 20, 2017
%  for use in Math 610 Lab
function [cell_matrix, area] = assembleLocalStiffness(localnodes, FE_at_Quad, Quad)
% Given an element defined by nodes in localnodes, we assemble the local stiffness
% term associated with this element.
%
% In 1D, the element is an interval, so there are p=2 nodes on the interval
% and p=2 linear basis functions corresponding to those nodes
%
% Input:
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
% cell_matrix = [pxp] matrix contributions from all interactions of basis function
% area        = [1x1] length of interval
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
% \int_{T} grad_lambda_i . grad_lambda_j dx
%   = \int_{That} (B'\grad lambdahat_i) . (B'\grad lamdbahat_j) |det(B)| dxhat
%


%affine map info
mat_B = (localnodes(2,:) - localnodes(1,:))';
det_B = abs(det(mat_B));
    
% precompute various things at the quadrature points if needed
% q_points = Quad.xhat*mat_B' + repmat(localnodes(1,:),Quad.nq,1); % [nqx2] list of quadrature points in cell 


% preallocate space for cell_matrix
cell_matrix = zeros(2,2);
area = det_B;

for q_index = 1:Quad.nq % run through quadrature points on element
    
    grad_phi_at_q_point=(mat_B')\[FE_at_Quad.hat_phix(q_index,:)]; % [dim x p] = [dimxdim][dimxp] p basis func gradients of length dim
  
    % exploiting matlab's speed with vector operations, we avoid doing
    % the for loops to construct the local matrix and rhs.  Below in
    % the commented section is the same code using for loops.
    grad_phi_ij_matrix = grad_phi_at_q_point'*grad_phi_at_q_point;  % [pxp] = [grad phi_i . grad_phi_j ]_{ij} using outer product
%     phi_ij_matrix = FE_at_Quad.hat_phi(q_index,:)'*FE_at_Quad.hat_phi(q_index,:);  %[pxp] = [phi_i phi_j]_{ij} using outer product

    cell_matrix = cell_matrix + grad_phi_ij_matrix' * Quad.what(q_index) * det_B;   


%     % The following shows the local assembling with the for loop.
%     for j = 1:2
%         for k = 1:2
%             cell_matrix(j,k) = cell_matrix(j,k) ...
%                      + dot(grad_phi_at_q_point(:,k),grad_phi_at_q_point(:,j))*Quad.what(q_index)*det_B;
%         end
%     end
    

end

