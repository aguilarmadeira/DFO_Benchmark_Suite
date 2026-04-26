function [f] = multi_gaussian(x);
%
% Purpose:
%
%    Function multi_gaussian is the function described in
%    Montaz et al. (2005).
%
%    dim = 2
%    Minimum global value = -1.29695
%    Global minima = [-0.01356 -0.01356]'
%    Local minima (five)
%    Search domain: -2 <= x <= 2 (Montaz et al. (2005))
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
A =[0.5 0 0 0.1
    1.2 1 0 0.5
    1 0 -0.5 0.5
    1 -0.5 0 0.5
    1.2 0 1 0.5];
f = -sum(A(:,1).*exp(-((x(1)-A(:,2)).^2+(x(2)-A(:,3)).^2)./A(:,4).^2),1);
%
% End of multi_gaussian.
