function [f] = poissonian(x);
%
% Purpose:
%
%    Function poissonian is the function described in
%    Brachetti et al. (1997).
%
%    dim = 2
%    Minimum global value = -95.28
%    Global minima
%    Local minima
%    Search domain: 1 <= x(1) <= 21 
%                   1 <= x(2) <= 8 (Brachetti et al. (1997))
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
n      = [5 2 4 2 7 2 4 5 4 4 15 10 8 15 5 6 3 4 5 2 6]';
i      = [1:21]';
lambda = 2+5.*exp(-0.5.*(((i-x(1))./x(2)).^2))+3;
f      = sum(-lambda + n.*log(lambda),1);
%
% End of poissonian.
