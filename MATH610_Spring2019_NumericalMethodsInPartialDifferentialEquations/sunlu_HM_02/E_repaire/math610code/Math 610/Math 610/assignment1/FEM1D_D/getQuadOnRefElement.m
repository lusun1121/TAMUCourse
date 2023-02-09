%  Written by Spencer Patty
%  Jan 20, 2017
%  for use in Math 610 Lab
function [Quad] = getQuadOnRefElement()
%
% return a structure containing the quadrature for a bulk 1D element
% (interval)
%
% To the ambitions student, you can add a parameter as input to choose
% between multiple quadratures that you have defined below.
%

% use a 2 pt quadrature exact for polynomials of degree 3
Quad.nq = 2;
% quadrature weights [nqx1]
Quad.what = [ 0.5; 0.5];
% quadrature points [nqxdim] = [nqx1]
Quad.xhat = [ 0.5 - 0.5/sqrt(3); 0.5+0.5/sqrt(3) ];

% % use a 1 pt quadrature exact for polynomials of degree 1
% Quad.nq = 1;
% % quadrature weights [nqx1]
% Quad.what = [ 1.0 ];
% % quadrature points [nqxdim] = [nqx1]
% Quad.xhat = [ 0.5];