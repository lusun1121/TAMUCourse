function [cell_matrix, local_rhs,area] = local_assemblyL2(localnodes,f,FE_at_Quad, Quad, p)

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

%
% x = F(xhat) 
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
mat_B = [(localnodes(2,:) - localnodes(1,:))',(localnodes(3,:) - localnodes(1,:))']; %[(dim)x(dim)]
det_B = det(mat_B);  

q_points = (mat_B*Quad.xhat' + repmat(localnodes(1,:)',1,Quad.nq))'; % [nqxdim] list of real quadrature points in cell not ref element
rhs_vals=zeros(Quad.nq,1);
for t=1:Quad.nq
rhs_vals(t) = f(q_points(t,:));
end

% preallocate space for cell_matrix
cell_matrix = zeros(3,3);
local_rhs = zeros(3,1);
area = abs(det_B);

for q_index = 1:Quad.nq % run through quadrature points on element
    
    
    phi_ij_matrix = FE_at_Quad.hat_phi(q_index,:)'*FE_at_Quad.hat_phi(q_index,:);  %3x3
    
    cell_matrix = cell_matrix + phi_ij_matrix...
                              * Quad.what(q_index) ...
                              area;%* abs(det_B);
                          
%   cell_matrix = cell_matrix + phi_ij_matrix...
%                             * Quad.what(q_index) ...
%                             * det_B;

    local_rhs = local_rhs + rhs_vals(q_index)...               % f(q_point)
                          * FE_at_Quad.hat_phi(q_index,:)' ... % bases on reference element at q_point [(p+1)x1]
                          * Quad.what(q_index)...              % quadrature weight
                          area;%* abs(det_B);                             % reference area-element transform   


end