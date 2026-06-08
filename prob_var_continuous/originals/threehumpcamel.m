function [f] = threehumpcamel(x);
%
% Purpose:
%
%    Function threehumpcamel is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = 0
%    Global minima = [0 0]'
%    Local minima (three)
%    Search domain: -5 <= x <= 5 (Montaz et al. (2005))
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
f = 2*x(1)^2-1.05*x(1)^4+1/6*x(1)^6+x(1)*x(2)+x(2)^2;
%
% End of threehumpcamel.
