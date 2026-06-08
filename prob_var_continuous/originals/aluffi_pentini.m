function [f] = aluffi_pentini(x);
%
% Purpose:
%
%    Function aluffi_pentini is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -0.3523
%    Global minima = [-1.0465 0]'
%    Local minima (two)
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
f = 0.25*x(1)^4-0.5*x(1)^2+0.1*x(1)+0.5*x(2)^2;
%
% End of aluffi_pentini.
