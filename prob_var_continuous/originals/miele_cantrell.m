function [f] = miele_cantrell(x);
%
% Purpose:
%
%    Function miele_cantrell is the function described in
%    Brachetti et al. (1997).
%
%    dim = 4
%    Minimum global value = 0
%    Global minima = [0 1 1 1]'
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
f = (exp(x(1))-x(2))^4 + 100*(x(2)-x(3))^6 + tan(x(3)-x(4))^4 + x(1)^8;
%
% End of miele_cantrell.
