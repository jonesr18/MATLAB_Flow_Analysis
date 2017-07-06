function out=lin2logicle(S)
% LIN2LOGICLE(S) converts a vector of values from a linear scale to a logicle scale 
%   OUT = LIN2LOGICLE(S) returns a vector of values corresponding to the values 
%   in X transformed from a linear scale to a logicle one using
%   the conversion function described by Parks, et al. "A New Logicle
%   Display Method Avoids Deceptive Effects of Logarithmic Scaling for Low
%   Signals and Compensated Data" Equation (5), where
%     S - a vector of linear 'raw' values
%     out - a vector of logicle converted values from S
%
%   Since Equation 5 cannot be explicitly solved for X, a spline is fitted to
%   X and S data using logicle2lin.m.  This spline is saved in the same folder 
%   and used to evaluated the values in S.
%
%   Example:
%       out = lin2logicle(linspace(0,1000));
%
%   Written by
%   Breanna Stillo
%   bstillo@mit.edu
%   Last Updated: 2014-10-14




% if exist('biexpfit.m','file')
%     load('biexpfit')
% else
% 
%     X=[linspace(0,4.5,1000)];
%     Y=logicle2lin(X);
% 
%     p=spline(Y,X);
% 
%     save('biexpfit','p')
% end
% 
% out=ppval(p,S);


T=262144;   % The top data value
M=4.5;      % Breadth of the display in decades
r=-150;    % Negative range reference value

W=(M-log10(T/abs(r)))/2;

A=0;

out=M.*logicleTransform(S,T,W,M,A);

end