function [z,mu,sigma] = zscore2(varargin)
%ZSCORE2 Standardized z score.
%   Z = ZSCORE2(X) returns a centered, scaled version of X, the same size as X.
%   For vector input X, Z is the vector of z-scores (X-MEAN(X)) ./ STD(X). For
%   matrix X, z-scores are computed using the mean and standard deviation
%   along each column of X.  For higher-dimensional arrays, z-scores are
%   computed using the mean and standard deviation along the first
%   non-singleton dimension.
%
%   The columns of Z have sample mean zero and sample standard deviation one
%   (unless a column of X is constant, in which case that column of Z is
%   constant at 0).
%
%   [Z,MU,SIGMA] = ZSCORE2(X) also returns MEAN(X) in MU and STD(X) in SIGMA.
%
%   [...] = ZSCORE2(X,1) normalizes X using STD(X,1), i.e., by computing the
%   standard deviation(s) using N rather than N-1, where N is the length of
%   the dimension along which ZSCORE works.  ZSCORE2(X,0) is the same as
%   ZSCORE2(X).
%
%   [...] = ZSCORE2(X,W,DIM) standardizes X by working along the dimension
%   DIM of X. Pass in W==0 to use the default normalization by N-1, or 1
%   to use N.
%	
%	ZSCORE2(...,NANFLAG) specifies how NaN (Not-A-Number) values are treated. 
%	The default is 'includenan':
%
%	'includenan' - the z-score of a vector containing NaN 
%                  values is also NaN.
%   'omitnan'    - elements of X or W containing NaN values are ignored.
%                  If all elements are NaN, the result is NaN.
%
%   See also ZSCORE, MEAN, STD.
%
%	Modified from ZSCORE by Ross Jones 2019-08-14 to have the NANFLAG option

%   Copyright 1993-2015 The MathWorks, Inc. 


% [] is a special case for std and mean, just handle it out here.
x = varargin{1};
if isequal(x,[]), z = x; return; end

if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
elseif isnumeric(varargin{3})
	dim = varargin{3};
	assert(dim <= ndims(x), 'dims argument higher than dimensionality of x');
end

if ischar(varargin{end})
	nanflag = varargin{end};
else
	nanflag = 'includenan';
end

% Compute X's mean and sd, and standardize it
mu = mean(x,dim,nanflag);
sigma = std(varargin{:});
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);