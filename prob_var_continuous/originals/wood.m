function [f] = wood(x);
%
% Purpose:
%
%    Function wood is the function described in
%    Brachetti et al. (1997).
%
%    dim = 4
%    Minimum global value = 0
%    Global minima = [1 1 1 1]'
%    Local minima 
%    Search domain: -10 <= x <= 10 (Brachetti et al. (1997))
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
f = 100*(x(1,1)^2 - x(2,1))^2 + (1-x(1,1))^2 + 90*(x(3,1)^2-x(4,1))^2 +...
   (1-x(3,1))^2 + 10.1*((1-x(2,1))^2 + (1-x(4,1))^2) + 19.8*(1-x(2,1))*(1-x(4,1));
%
% End of wood.
