function [f] = becker_lago(x);
%
% Purpose:
%
%    Function becker_lago is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0
%    Global minima = [-5 -5]', [-5 5]', [5 -5]', [5 5]'
%    Local minima 
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
f = (abs(x(1))-5)^2+(abs(x(2))-5)^2;
%
% End of becker_lago.
