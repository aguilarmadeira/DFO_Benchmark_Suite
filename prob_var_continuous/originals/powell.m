function [f] = powell(x);
%
% Purpose:
%
%    Function powell is the function described in
%    Montaz et al. (2005).
%
%    dim = 4
%    Minimum global value = 0
%    Global minima = zeros(4,1)
%    Local minima (equal to global minima)
%    Search domain: -10 <= x <= 10 (Montaz et al. (2005))
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
f = (x(1,1) + 10*x(2,1))^2 + 5*(x(3,1) - x(4,1))^2 + (x(2,1) -...
    2*x(3,1))^4 + 10*(x(1,1) - x(4,1))^4;
%
% End of powell.
