function [f] = sinusoidal(x);
%
% Purpose:
%
%    Function sinusoidal is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = -(A+1)
%    Global minima = (90+z)*ones(n,1)
%    Local minima (several)
%    Search domain: 0 <= x <= 180 (Montaz et al. (2005))
%    Cases considered: n = 10, 20
%                      A = 2.5;B = 5;z = 30;Montaz et al. (2005)
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
A = 2.5;
B = 5;
z = 30;
f = -(A*prod(sind(x-z),1)+prod(sind(B.*(x-z)),1));
%
% End of sinusoidal.
