function [f] = hosaki(x);
%
% Purpose:
%
%    Function hosaki is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -2.3458
%    Global minima = [4 2]'
%    Local minima (two)
%    Search domain: 0 <= x(1) <= 5
%                   0 <= x(2) <= 6 (Montaz et al. (2005))
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
f = (1-8*x(1)+7*x(1)^2-7/3*x(1)^3+1/4*x(1)^4)*x(2)^2*exp(-x(2));
%
% End of hosaki.
