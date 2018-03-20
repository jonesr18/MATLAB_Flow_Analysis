function r = log10nrnd(mu,sigma,varargin)
%LOG10NRND Random arrays from the lognormal distribution.
%   R = LOG10NRND(MU,SIGMA) returns an array of random numbers generated from 
%   the lognormal distribution with parameters MU and SIGMA.  MU and SIGMA 
%   are the mean and standard deviation, respectively, of the associated 
%   normal distribution.  The size of R is the common size of MU and SIGMA 
%   if both are arrays.  If either parameter is a scalar, the size of R is 
%   the size of the other parameter.
%
%   R = LOG10NRND(MU,SIGMA,M,N,...) or R = LOG10NRND(MU,SIGMA,[M,N,...])
%   returns an M-by-N-by-... array.
%
%	This function uses log10 transformations rather than log (as in LOGNPDF). 
%	
%   The mean and variance of a lognormal random variable with parameters MU
%   and SIGMA are
%
%      M = 10.^(MU + SIGMA^2/2)
%      V = 10.^(2*MU + SIGMA^2) * (10.^(SIGMA^2) - 1)
%
%   Therefore, to generate data from a lognormal distribution with mean M and
%   Variance V, use
%
%      MU = log10(M^2 / sqrt(V+M^2))
%      SIGMA = sqrt(log10(V/M^2 + 1))
%
%   See also LOGNRND, LOG10NPDF, LOG10NSTAT, 

%   LOG10NRND uses a transformation of a normal random variable.

%   References:
%      [1] Marsaglia, G. and Tsang, W.W. (1984) "A fast, easily implemented
%          method for sampling from decreasing or symmetric unimodal density
%          functions", SIAM J. Sci. Statist. Computing, 5:349-359.
%      [2] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2015 The MathWorks, Inc.
%	--> Copied from lognrnd and edited by Ross Jones


if nargin < 2
    error(message('stats:lognrnd:TooFewInputs'));
end

[err, sizeOut] = statsizechk(2,mu,sigma,varargin{:});
if err > 0
    error(message('stats:lognrnd:InputSizeMismatch'));
end

% Return NaN for elements corresponding to illegal parameter values.
sigma(sigma < 0) = NaN;

ty = internal.stats.dominantType(mu, sigma);
r = 10.^(randn(sizeOut,'like',ty) .* sigma + mu);
