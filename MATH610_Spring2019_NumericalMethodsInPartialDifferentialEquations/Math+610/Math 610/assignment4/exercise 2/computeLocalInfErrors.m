function [local_Linf_error] = computeLocalInfErrors(localnodes, uh_at_localnodes, u_exact, Quad, FE_at_Quad,p)
% Given an element defined by nodes in localnodes, we compute the Linfinity
% error compared to the exact solutions on the local element.
%
% In 1D, the element is an interval, so there are (dim+1)=2 nodes on the interval
% and p+1 lagrange basis functions of degree p
%
% Input:
% localnodes: [(dim+1) x dim] is a column vector of nodes on element (each node is a "row" vector)
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

mat_B = [(localnodes(2,:) - localnodes(1,:));(localnodes(3,:) - localnodes(1,:))]; %[(dim)x(dim)]

q_points = Quad.xhat*mat_B + repmat(localnodes(1,:),Quad.nq,1); % [nqx(dim+1)] list of quadrature points in cell 

% evaluate solutions and gradients at quadrature points
u_exact_at_q_points=zeros(Quad.nq,1);
for t=1:Quad.nq
u_exact_at_q_points(t) = u_exact(q_points(t,:)); %[nqx1]
end
uh_at_q_points = FE_at_Quad.hat_phi*uh_at_localnodes;  %[nqx1] = [nqx(p+1)] [(p+1)x1]

% assemble error analysis
local_Linf_error = 0;

for q_index = 1:Quad.nq % run through quadrature points on element
    diff_val = uh_at_q_points(q_index) - u_exact_at_q_points(q_index); %[1x1] 
    local_Linf_error = max(local_Linf_error, abs(diff_val));
end




