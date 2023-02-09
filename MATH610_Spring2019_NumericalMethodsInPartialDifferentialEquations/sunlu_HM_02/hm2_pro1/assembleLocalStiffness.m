%  Written by Spencer Patty
%  Jan 20, 2017
%  for use in Math 610 Lab
function [cell_matrix, area] = assembleLocalStiffness(localnodes, FE_at_Quad, Quad, p)
% Given an element defined by nodes in localnodes, we assemble the local stiffness
% term associated with this element.
%
% In 1D, the element is an interval, so there are (dim+1)=2 nodes on the interval
% and p+1 lagrange basis functions of order p on each element
%
% Input:
% localnodes: [(dim+1)xdim] is a column vector of vertices on element (each vertex is a "row" vector)
% Quad.nq: [1x1] number of quadrature points on reference element
% Quad.what: [nqx1] weights on reference element for each quad point
% Quad.xhat: [nqxdim] coordinates of quad point on reference element
% FE_at_Quad.hat_phi: [nqx(p+1)] p+1 basis functions on reference element
%                                evaluated at the Quad.nq quad points
% FE_at_Quad.hat_phix: [nqx(p+1)] x derivative of p+1 basis functions on 
%                                reference element evaluated at the Quad.nq 
%                                quad points
%
% Output:
% cell_matrix = [p+1xp+1] matrix contributions from all interactions of basis function
% area        = [1x1] length of interval
%
% Let localnodes = [v_1; ... v_{dim+1}] where each v_i [1xdim] is the ith node
% on local cell T
%
% In 1D, localnodes = [v_1; v_2] = [2x1]
% In 2D, localnodes = [v_1; v_2; v_3]  = [3x2], etc

% 1D Conversion from reference element, That = (0,1) to T = (v1, v2)
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
% \int_{T} grad_lambda_i . grad_lambda_j dx
%   = \int_{That} (B\grad lambdahat_i) . (B\grad lamdbahat_j) |det(B)| dxhat
%


%affine map info
%mat_B = (localnodes(2,:) - localnodes(1,:))';
mat_B = (localnodes(3,:) - localnodes(1,:))';
det_B = abs(det(mat_B));
inv_B = 1/mat_B;

% precompute various things at the quadrature points if needed
% q_points = Quad.xhat*mat_B + repmat(localnodes(1,:),Quad.nq,1); % [nqxdim] list of quadrature points in cell 


% preallocate space for cell_matrix
cell_matrix = zeros(p+1,p+1);
%cell_matrix.fi0 = zeros(p+1,p+1);
area = det_B;

for q_index = 1:Quad.nq % run through quadrature points on element
    
    grad_phi_at_q_point = inv_B*FE_at_Quad.hat_phix(q_index,:); % [dimx(p+1)] = [dimxdim][dimx(p+1)] p+1 basis func gradients of length dim
  
    % assume j index is for solution basis functions and i index is for 
    % test functions basis.  We should return A_{ij} to get the indexing 
    % correct for non symmetric problems.
    
    % exploiting matlab's speed with vector operations, we avoid doing
    % the for loops to construct the local matrix and rhs.  Below in
    % the commented section is the same code using for loops.
    grad_phi_ij_matrix = grad_phi_at_q_point'*grad_phi_at_q_point;  % [(p+1)x(p+1)] = [(p+1)xdim][dimx(p+1)] = [grad phi_i . grad_phi_j ]_{ij} using outer product
    %phi_ij_matrix = FE_at_Quad.hat_phi(q_index,:)'*FE_at_Quad.hat_phi(q_index,:);  %[(p+1)x(p+1)] = [(p+1)x1][1x(p+1)] = [phi_i phi_j]_{ij} using outer product

    cell_matrix = cell_matrix + grad_phi_ij_matrix ...
                              * Quad.what(q_index) ...
                              * det_B; 
   
    %cell_matrix.fi0 = cell_matrix.fi0 + phi_ij_matrix ...
     %                         * Quad.what(q_index) ...
     %                         * det_B; 


%     % The following shows the local assembling with the for loop.
%     for i = 1:(p+1) // test functions
%         for j = 1:(p+1)
%             cell_matrix(i,j) = cell_matrix(i,j) ...
%                      + dot(grad_phi_at_q_point(:,i),grad_phi_at_q_point(:,j))*Quad.what(q_index)*det_B;
%         end
%     end
    

end

