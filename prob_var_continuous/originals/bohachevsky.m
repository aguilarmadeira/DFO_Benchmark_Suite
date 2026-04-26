function [f] = bohachevsky(x);
%
% Purpose:
%
%    Function bohachevsky is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0
%    Global minima = [0 0]'
%    Local minima 
%    Search domain: -50 <= x <= 50 (Montaz et al. (2005))
%
% Input:  
%
%         x (point given by the optimizer).
%
% Output: 
%
%         f (function value at x).
%
% Written by A. L. Custodio and J. F. A. Madeira.
%
% Version June 2012.
%
%
f = x(1)^2+2*x(2)^2-0.3*cos(3*pi*x(1))*cos(4*pi*x(2))+0.3;
%
% End of bohachevsky.
