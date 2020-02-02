function gsd = geostd(x, flag, dim, nanflag)
% For vectors, GEOSTD(X) calculates the geometric standard deviation of
% input X.  For  matrices GEOSTD calculates the geometric standard
% deviation over each column of X.  For N-D arrays, GEOSTD calculates the
% geometric standard deviation over the first non-singleton dimension of X.
% 
% 
% GEOSTD(X,1) normalizes by N and produces the square root of the second
% moment of the sample about its mean. GEOSTD(X,0) is the same as GEOSTD(X).
% 
% GEOSTD(X, [], DIM) calculates the geometric standard deviation along
% dimension DIM of X.  Use a flag of 0 to normalise by factor (N-1), or a
% flag of 1 to normalise by factor N.
% 
% GEOSTD(X, [], DIM, NANFLAG) specifies how NaN (Not-A-Number) values are treated.
% The default is 'includenan':
% 
% 'includenan' - the standard deviation of a vector containing NaN 
%                values is also NaN.
% 'omitnan'    - elements of X or W containing NaN values are ignored.
%                If all elements are NaN, the result is NaN.
% 
% NOTE: Class type checking and error handling are conducted within EXP,
% STD and LOG.
%
% 
% EXAMPLE:  X = 10*rand(5);
%                      geostd(X)
%                      ans =
%                     
%                         1.1858    1.8815    1.8029    4.1804    2.5704
% 
% 
%   Class support for input X:
%      float: double, single
% 
% 
%   See also GEOMEAN (stats toolbox), STD.
% 
% 
% $ Author: Richie Cotton $     $ Date: 2006/03/17 $
%
% Update Log:
%	2017-09-28		Ross Jones
%		Empty input now returns NaN output
%	2019-05-03		Ross Jones
%		Added option to omit NaNs from analysis
%		Changed default behavior to 'normal' for std normalization


% Basic error checking of inputs
if nargin < 1
    error('geostd:notEnoughInputs', 'This function requires at least one input.');
elseif any(x(:) < 0)
    error(geostd:badData', 'All data values must be positive.');
end

% Return NaN for empty inputs
if isempty(x)
	gsd = NaN;
	return;
end

% Setup default flag where required
if nargin < 2 || isempty(flag)
    flag = 0;
end

% If dimension is not specified, find first non-singleton dimension
if nargin < 3 || isempty(dim)
    dim = find(size(x) ~= 1, 1);
    if isempty(dim)
        dim = 1;
    end
end

% If nanflag is not specified, include NaNs
if nargin < 4 || isempty(nanflag)
	nanflag = 'includenan';
end

% Turn off warnings regarding log of zero, since this is an artifact of the
% technique use to calculate the gsd
lozwarning = warning('off', 'MATLAB:log:logOfZero');

% Calculate geometric std dev using 
% "log of geometric std dev of data = arithmetic std dev of log of data"
gsd = exp(std(log(x), flag, dim, nanflag));

% Reset warning value
warning(lozwarning);