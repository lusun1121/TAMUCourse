%  Written by Spencer Patty
%  Jan 20, 2017
%  for use in Math 610 Lab
function [FE_at_quad] = feEval( Quad, p )
% Return shape value and shape gradients on quadrature points
% in the reference element
%
%
% Outputs the p+1 Lagrange basis functions of order p on ref element [0,1]
% evaluated at the quadrature points along with the derivative components
% of each basis function to construct the gradients.
%
% input:
%   Quad.nq = number of quadrature points
%   Quad.xhat: [nqxdim] list of nq quadrature points [x]
%
% output: 
%   hat_phi: [nqx(p+1)] shape value 
%   hat_phix: [nqx(p+1)] shape gradient x component
%

% order of basis functions
% 1) left dof
% 2) right dof
% 3) interior dofs from left to right

x = Quad.xhat(:,1);
% y = Quad.xhat(:,2);

if (p == 1)
    %                   phi_0, phi_1
    % Linear elements :   1-x,  x
    % x derivative    :   -1,   1
    FE_at_quad.hat_phi  = [1-x, x];
    FE_at_quad.hat_phix = [-1*ones(Quad.nq,1), 1*ones(Quad.nq,1)];
elseif (p==2)
   % construct 3 basis functions on [0,1]:  phi_0,  phi_1,  phi_{1/2}
   
   % construct 3 x-derivatives of basis functions on [0,1]:  dphi_0/dx, dphi_1/dx,  dphi_{1/2}/dx
   FE_at_quad.hat_phi = [(2*x-1).*(x-1), 4*x.*(1-x),(2*x-1).*x];
   FE_at_quad.hat_phix = [4*x-3,   -8*x+4, 4*x-1];   
elseif (p==3)
   % construct 4 basis functions on [0,1]:  phi_0,  phi_1,  phi_{1/3}, phi_{2/3}
   
   % construct 4 x-derivatives of basis functions on [0,1]:  dphi_0/dx, dphi_1/dx,  dphi_{1/3}/dx , dphi_{2/3}/dx
    
end



