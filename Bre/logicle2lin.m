function out=logicle2lin(X)
% LOGICLE2LIN(X) converts a vector of values from a logicle scale to a linear scale 
%   OUT = LOGICLE2LIN(X) returns a vector of values corresponding to the values 
%   in X transformed from a logicle scale to a linear one using
%   the conversion function described by Parks, et al. "A New Logicle
%   Display Method Avoids Deceptive Effects of Logarithmic Scaling for Low
%   Signals and Compensated Data" Equation (5), where
%     X - a vector of logicle values
%     out - a vector of linear converted values from X
%
%   Example:
%       out = logicle2lin(linspace(0,4));
%
%   Written by
%   Breanna Stillo
%   bstillo@mit.edu
%   Last Updated: 2014-10-14


T=262144;   % The top data value
M=4.5;      % Breadth of the display in decades
r=-150;    % Negative range reference value

W=(M-log10(T/abs(r)))/2;
% 
% syms p
% P=solve(W==2*p*log10(p)/(p+1),p);
% 
% posI=X>=W;
% S=@(X) subs(T.*10^-(M-W).*(10.^(X-W)-P^2.*10.^(-(X-W)./P)+P^2-1));
% 
% out=eval(S(X).*posI-S(2*W-X).*(1-posI));



A=0;
out=logicleInverseTransform(X./M,T,W,M,A);

end