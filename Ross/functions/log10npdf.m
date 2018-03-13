function y = log10npdf(x, mu, sigma, param)
%LOG10NPDF Lognormal probability density function (pdf).
%   Y = LOG10NPDF(X,MU,SIGMA,PARAM) returns values at X of the lognormal pdf 
%   with distribution parameters MU and SIGMA. MU and SIGMA are the mean and 
%   standard deviation, respectively, of the associated normal distribution.  
%   The size of Y is the common size of the input arguments. A scalar input 
%   functions as a constant matrix of the same size as the other inputs.
%
%	PARAM indicates how to interpret MU. Defualt is 'mu', but can accept
%	'mean', 'median', and 'mode', which will automatically adjust the value 
%	of MU to one of these values:
%		mean	= exp(MU + SIGMA^2 / 2)
%		median	= exp(MU)
%		mode	= exp(MU - SIGMA^2)
%	
%	This function uses log10 transformations rather than log (as in LOGNPDF). 
%	
%   Default values for MU and SIGMA are 0 and 1 respectively.
%
%   See also LOGNPDF, LOG10NRND, LOG10NSTAT.

%   References:
%      [1] Evans, M., Hastings, N., and Peacock, B. (1993) Statistical
%          Distributions, 2nd ed., Wiley, 170pp.

%   Copyright 1993-2006 The MathWorks, Inc.
%	--> Copied from lognpdf and edited by Ross Jones


if nargin<1
    error(message('stats:lognpdf:TooFewInputs'));
end
if nargin < 2
    mu = 0;
end
if nargin < 3
    sigma = 1;
end
if nargin < 4
	param = 'mean';
else
	validatestring(param, {'mean', 'median', 'mode', 'mu'}, mfilename, 'param', 4);
end

% Return NaN for out of range parameters.
sigma(sigma <= 0) = NaN;

% Negative data would create complex values, potentially creating spurious
% NaNi's in other elements of y.  Map them, and zeros, to the far right
% tail, whose pdf will be zero.
x(x <= 0) = Inf;

switch param
	case {'mu', 'median'}
		% Do nothing: mu = mu
	case 'mean'
		mu = mu + sigma^2 / 2;
	case 'mode'
		mu = mu - sigma^2;
end

try
    y = 10.^(-0.5 * ((log10(x) - mu)./sigma).^2) ./ (x .* sqrt(2*pi) .* sigma);
catch
    error(message('stats:lognpdf:InputSizeMismatch'));
end