function y = log10poisspdf(x, lambda)
%LOG10CPOISSPDF Base-10 Log Poisson probability density function.
%   Y = LOG10POISSPDF(X,LAMBDA) returns the base-10 log Poisson 
%   probability density function with parameter LAMBDA at the values in X.
%
%   The size of Y is the common size of X and LAMBDA. A scalar input   
%   functions as a constant matrix of the same size as the other input.    
%
%   See also POISSPDF, POISSRND, POISSTAT

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

k = find(x >= 0 & x == round(x) & lambda > 0);

if ~isempty(k)
    x = x(k);
    lambda = lambda(k);
 
    smallx = x<=lambda*realmin;
    y(k(smallx)) = exp(-lambda(smallx));
    
    largex = lambda<x*realmin;
    y(k(largex)) = exp(-lambda(largex) + log10(x(largex)).*log(lambda(largex)) ...
        - gammaln(log10(x(largex))+1) - log(x(largex)));
    
    other = ~smallx & ~largex;
%     lnsr2pi = 0.9189385332046727; % log(sqrt(2*pi))
%     y(k(other)) = exp(-lnsr2pi -0.5*log(log10(x(other))) - stirlerr(log10(x(other))) ...
%         - binodeviance(log10(x(other)),lambda(other)) - x); % log10(x) in binodeviance?
	y(k(other)) = (lambda(other).^log10(x(other)) .* exp(-lambda(other))) ./ ...
			gamma(log10(x(other)) + 1) ./ x(other);
end
