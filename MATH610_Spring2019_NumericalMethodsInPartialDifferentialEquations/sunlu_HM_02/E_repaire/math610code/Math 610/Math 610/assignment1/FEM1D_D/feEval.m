%  Written by Spencer Patty
%  Jan 20, 2017
%  for use in Math 610 Lab
function [FE_at_quad] = feEval( Quad )
% Return shape value and shape gradients on quadrature points
% in the reference element
% The following is the linear case. 
% input:
%   Quad.nq = number of quadrature points
%   Quad.xhat: [nqxdim] list of nq quadrature points [x]
%
% output:  p = number of basis functions
%   hat_phi: [nqxp] shape value 
%   hat_phix: [nqxp] shape gradient x component
%

% Linear elements :  1-x,  x
% x derivative    :   -1,  1


FE_at_quad.hat_phi  = [1-Quad.xhat(:,1), Quad.xhat(:,1)];
FE_at_quad.hat_phix = [-1*ones(Quad.nq,1), 1*ones(Quad.nq,1)];


