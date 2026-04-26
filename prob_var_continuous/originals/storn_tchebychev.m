function [f] = storn_tchebychev(x);
%
% Purpose:
%
%    Function storn_tchebychev is the function described in
%    Montaz et al. (2005).
%
%    dim = n
%    Minimum global value = 0
%    Global minima: [128 0 -256 0 160 0 -32 0 1]', for n = 9
%    Local minima
%    Search domain: -128 <= x <= 128, for n = 9
%                   -32768 <= x <= 32768, for n = 17 (Montaz et al. (2005))
%    Cases considered: n = 9, d = 72.661, m = 60
%                      n = 17, d = 10558.145, m = 100 Montaz et al. (2005)
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
dim = size(x,1);
if dim == 9
    d = 72.6066;
    m = 60;
else
    if dim == 17
       d = 10558.145;
       m = 100;
    else
        return
    end
end
auxu = (dim-1-[0:dim-1])';
u    = sum(1.2.^auxu.*x,1);
v    = sum((-1.2).^auxu.*x,1);
if u<d
    p1 = (u-d)^2;
else
    p1 = 0;
end
if v<d
    p2 = (v-d)^2;
else
    p2 = 0;
end
w = sum((repmat(2.*[0:m]./m-1,dim,1).^repmat(auxu,1,m+1)).*repmat(x,1,m+1),1);
for i = 1:m
    if w(i) >1
        paux(i,1) = (w(i)-1)^2;
    else
        if w(i)<-1
            paux(i,1) = (w(i)+1)^2;
        else
            paux(i,1) = 0;
        end
    end
end
f = p1+p2+sum(paux,1);
%
% End of storn_tchebychev.
