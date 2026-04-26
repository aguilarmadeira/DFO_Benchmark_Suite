function [f] = meyer_roth(x);
%
% Purpose:
%
%    Function meyer_roth is the function described in
%    Brachetti et al. (1997).
%
%    dim = 3
%    Minimum global value = 0.4 * 10^-4
%    Global minima = [3.13 15.16 0.78]'
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
T = [1 2 1 2 0.1]';
V = [1 1 2 2 0]';
Y = [0.126 0.219 0.076 0.126 0.186]';
f = sum(((x(1)*x(3).*T)./(1+x(1).*T+x(2).*V)-Y).^2,1);
%
% End of meyer_roth.
