function y = cpoisspdf(x, lambda)
%CPOISSPDF Continuous Poisson probability density function.
%   Y = CPOISSPDF(X,LAMBDA) returns the continuous Poisson probability 
%   density function with parameter LAMBDA at the values in X.
%
%   The size of Y is the common size of X and LAMBDA. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   See also POISSPDF, CPOISSRND, POISSTAT

%   References:
%      [1]  M. Abramowitz and I. A. Stegun, "Handbook of Mathematical
%      Functions", Government Printing Office, 1964, 26.1.22. 
%      [2]  C. Loader, "Fast and Accurate Calculations of Binomial
%      Probabilities", 2000.

%   Copyright 1993-2015 The MathWorks, Inc.
%	--> Copied from poisspdf and edited by Ross Jones


if nargin <  2, 
    error(message('stats:poisspdf:TooFewInputs')); 
end

[errorcode, x, lambda] = distchck(2,x,lambda);

if errorcode > 0
    error(message('stats:poisspdf:InputSizeMismatch'));
end

y = zeros(size(x),'like',internal.stats.dominantType(x,lambda)); % Single if x or lambda is

if ~isfloat(x)
   x = double(x);
end

y(lambda < 0) = NaN;
y(isnan(x) | isnan(lambda)) = NaN;
y(x==0 & lambda==0) = 1;

k = find(x >= 0 & lambda > 0);

if ~isempty(k)
    x = x(k);
    lambda = lambda(k);
 
    smallx = x<=lambda*realmin;
    y(k(smallx)) = exp(-lambda(smallx));
    
    largex = lambda<x*realmin;
    y(k(largex)) = exp(-lambda(largex) + x(largex).*log(lambda(largex)) ...
        - gammaln(x(largex)+1));
    
    other = ~smallx & ~largex;
    lnsr2pi = 0.9189385332046727; % log(sqrt(2*pi))
    y(k(other)) = exp(-lnsr2pi -0.5*log(x(other)) - stirlerr(x(other)) ...
        - binodeviance(x(other),lambda(other)));
end
