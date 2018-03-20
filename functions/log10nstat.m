function [m,v]= log10nstat(mu,sigma)
%LOG10NSTAT Mean and variance for the lognormal distribution.
%   [M,V] = LOG10NSTAT(MU,SIGMA) returns the mean of and variance of the 
%   lognormal distribution with parameters MU and SIGMA.  MU and SIGMA are 
%   the mean and standard deviation, respectively, of the associated normal 
%   distribution.  The sizes of M and V are the common size of the input 
%   arguments.  A scalar input functions as a constant matrix of the same 
%   size as the other inputs.
%
%	This function uses log10 transformations rather than log (as in LOGNPDF). 
%	
%   See also LOGNSTAT, LOG10NPDF, LOG10NRND.

%   References:
%      [1] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2006 The MathWorks, Inc.
%	--> Copied from lognstat and edited by Ross Jones


if nargin < 2
    error(message('stats:lognstat:TooFewInputs'));
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

s2 = sigma .^ 2;

%   The inverse transformation is
%
%      MU = log10(M^2 / sqrt(V+M^2))
%      SIGMA = sqrt(log10(V/M^2 + 1))
try
    m = 10.^(mu + 0.5 * s2);
    v = 10.^(2*mu + s2) .* (10.^(s2)-1);
catch
    error(message('stats:lognstat:InputSizeMismatch'));
end