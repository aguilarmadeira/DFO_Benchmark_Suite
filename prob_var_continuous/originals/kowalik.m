function [f] = kowalik(x);
%
% Purpose:
%
%    Function kowalik is the function described in
%    Montaz et al. (2005).
%
%    dim = 4
%    Minimum global value = 3.0748*10^-4
%    Global minima = [0.192 0.190 0.123 0.135]'
%    Local minima 
%    Search domain: 0 <= x <= 0.42 (Montaz et al. (2005))
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
a = [0.1957 0.1947 0.1735 0.16 0.0844 0.0627 0.0456 0.0342 0.0323 0.0235 0.0246];
b = [0.25 0.5 1 2 4 6 8 10 12 14 16];
f = sum((a-x(1).*(1+x(2).*b)./(1+x(3).*b+x(4).*b.^2)).^2,2);
%
% End of kowalik.
