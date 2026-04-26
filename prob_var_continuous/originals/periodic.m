function [f] = periodic(x);
%
% Purpose:
%
%    Function periodic is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0.9
%    Global minima = [0 0]'
%    Local minima (49 with value 1)
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
f = 1+sin(x(1))^2+sin(x(2))^2-0.1*exp(-x(1)^2-x(2)^2);
%
% End of periodic.
