function [f] = cauchy(x);
%
% Purpose:
%
%    Function cauchy is the function described in
%    Brachetti et al. (1997).
%
%    dim = n
%    Minimum global value
%    Global minima
%    Local minima
%    Search domain: y(1) <= x <= y(dim) (Brachetti et al. (1997))
%    Cases considered: n = 4 Brachetti et al. (1997)
%                      n = 10 Brachetti et al. (1997)
%                      n = 25 Brachetti et al. (1997)
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
n = size(x,1);
if n == 4
    y = [3 7 12 17]';
else
    if n == 10
        y = [2 5 7 8 11 15 17 21 23 26]';
    else
        if n == 25
            y = [4.1 7.7 17.5 31.4 32.7 92.4 115.3 118.3 119 129.6 198.6];
            y = [y, 200.7 242.5 255 274.7 274.7 303.8 334.1 430 489.1 703.4 ];
            y = [y, 978 1656 1697.8 2745.6]';
        end
    end
end
f = -sum(log(pi)+log(1+(y-x).^2),1);
%
% End of cauchy.
