function [f] = rosenbrock(x);
%
% Purpose:
%
%    Function rosenbrock is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value = 0
%    Global minima = ones(n,1)
%    Local minima (equal to global minima)
%    Search domain: -1000 <= x <= 1000 (Brachetti et al. (1997))
%                   -2.048 <= x <= 2.048 (Storn and Price (1997))
%                   -5.12 <= x <= 5.12 (Huyer and Neumaier (2008))
%    Cases considered: n = 2 Huyer and Neumaier (2008)
%                            Brachetti et al. (1997)
%                      n = 10 Garcia-Palomares  et al. (2006)
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
n    = size(x,1);
xaux = x(2:n);
f    = sum((x(1:n-1) - 1).^2 + 100.*(x(1:n-1).^2 - xaux).^2,1);
%
% End of rosenbrock.
