function [f] = epistatic_michalewicz(x);
%
% Purpose:
%
%    Function epistatic_michalewicz is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value: -4.687 if n = 5
%                          -9.66 if n = 10
%    Global minima
%    Local minima (several)
%    Search domain: 0 <= x <= pi (Huyer and Neumaier (1999))
%    Cases considered: n = 5 Huyer and Neumaier (1999), ICEO
%                      n = 10 Huyer and Neumaier (1999), ICEO
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
even = 0;
for i=1:n-1
    if even == 0
       y(i,1) = x(i)*cos(pi/6)-x(i+1)*sin(pi/6);
       even = 1;
    else
       y(i,1) = x(i)*sin(pi/6)+x(i+1)*cos(pi/6);
       even = 0;
    end
end
y = [y;x(n)];
f = -sum(sin(y).*(sin([1:n]'.*(y.^2)./pi)).^20,1);
%
% End of epistatic_michalewicz.
