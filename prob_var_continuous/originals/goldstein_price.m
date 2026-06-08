function [f] = goldstein_price(x);
%
% Purpose:
%
%    Function goldstein_price is the function described in
%    Brachetti et al. (1997).
%
%    dim = 2
%    Minimum global value = 3
%    Global minima = [0 -1]'
%    Local minima (four)
%    Search domain: -2 <= x <= 2 (Brachetti et al. (1997))
%                                (Huyer and Neumaier (1999))
%                                (Huyer and Neumaier (2008))
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
f = 1+(x(1)+x(2)+1)^2*(19-14*x(1)+3*x(1)^2-14*x(2)+6*x(1)*x(2)+3*x(2)^2);
f = f * (30+(2*x(1)-3*x(2))^2*(18-32*x(1)+12*x(1)^2+48*x(2)-36*x(1)*x(2)+27*x(2)^2));
%
% End of goldstein_price.
