function [f] = dekkers_aarts(x);
%
% Purpose:
%
%    Function dekkers_aarts is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -24771
%    Global minima = [0 15]', [0 -15]'
%    Local minima (three)
%    Search domain: -20 <= x <= 20 (Montaz et al. (2005))
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
f = 10^5*x(1)^2+x(2)^2-(x(1)^2+x(2)^2)^2+1/10^5*(x(1)^2+x(2)^2)^4;
%
% End of dekkers_aarts.
