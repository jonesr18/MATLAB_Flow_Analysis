function [p, rnorm, r, exitflag, output, lambda, J] = lsqmultifit(x_cell, y_cell, mdl_cell, p0, lb, ub, varargin)
%LSQMULTIFIT least-squares curve fit of multiple data sets
%
%	A wrapper function for LSQCURVEFIT which allows simulatenous fitting of
%	multiple data sets with shared fitting parameters. See example below.
%
%	INPUT:
% 		x_cell,y_cell: Cell arrays containing the x,y vectors of the fitted
%					   data sets.
% 		mdl_cell: Cell array containing model functions for each data set.
% 		p0: Vector containing initial guess of the fitted parameters.
% 		options: Structure containing control parameters for LSQCURVEFIT (see
%				 help file on LSQCURVEFIT for more details).
%		Name,Value pairs: additional options passed to LSQCURVEFIT (see help
%						  file on LSQCURVEFIT for more details). 
%
%	OUTPUT:
%		p, rnorm, r, exitflag, output, lambda, J: Direct output from LSQCURVEFIT.
%
%	EXAMPLE:
% 		% Generate X vectors for both data sets
% 		x1 = 0:0.1:10;
% 		x2 = 0:1:10;
%		
% 		% Generate Y data with some noise
% 		y1 = cos(2*pi*0.5*x1).*exp(-x1/5) + 0.05*randn(size(x1));
% 		y2 = 0.5 + 2*exp(-(x2/5)) + 0.05*randn(size(x2));
% 
% 		% Define fitting functions and parameters, with identical
%		% exponential decay for both data sets
% 		mdl1 = @(p, x) cos(2 * pi * p(1) * x) .* exp(-x / p(2));
% 		mdl2 = @(p, x) p(4) + p(3) * exp(-( x / p(2)));
% 
% 		% Prepare input for LSQMULTIFIT and perform fitting
% 		x_cell = {x1, x2};
% 		y_cell = {y1, y2};
% 		mdl_cell = {mdl1, mdl2};
% 		p0 = [1, 1, 1, 1];
% 		[p, rnorm, r, exitflag, output, lambda, J] = ...
%					nlinmultifit(x_cell, y_cell, mdl_cell, p0);
%		
% 		% Calculate parameter confidence intervals
% 		ci = nlparci(p, r, 'Jacobian', J);
%
%	AUTHOR:
%		Ross Jones
%		jonesr18@mit.edu
%
%		Adapted from nlinmultifit:
%			Chen Avinadav
%			mygiga (at) gmail
%
	
	if ~exist('lb', 'var'), lb = -inf * ones(size(p0)); end
	if ~exist('ub', 'var'), ub = inf * ones(size(p0)); end
	
	if (~iscell(x_cell) && ~iscell(y_cell) && ~iscell(mdl_cell))
		[p, rnorm, r, exitflag, output, lambda, J] = lsqcurvefit(mdl_cell, p0, x_cell, y_cell, lb, ub, varargin{1:end});
		return;
	end
	
	num_curves = length(x_cell);
	if length(y_cell) ~= num_curves || length(mdl_cell) ~= num_curves
		error('Invalid input to LSQMULTIFIT');
	end
	
	x_vec = [];
	y_vec = [];
	mdl_vec = '@(p, x) [';
	mdl_ind1 = 1;
	mdl_ind2 = 0;
		
	for ii = 1:num_curves
		if length(x_cell{ii}) ~= length(y_cell{ii})
			error('Invalid input to LSQMULTIFIT');
		end
		
		x_vec = [x_vec; x_cell{ii}]; %#ok<AGROW>
		y_vec = [y_vec; y_cell{ii}]; %#ok<AGROW>
		mdl_ind2 = mdl_ind2 + size(x_cell{ii}, 1);
		mdl_vec = [mdl_vec, sprintf('mdl_cell{%d}(p, x(%d:%d, :)); ', ii, mdl_ind1, mdl_ind2)]; %#ok<AGROW>
		mdl_ind1 = mdl_ind1 + size(x_cell{ii}, 1);
	end
	mdl_vec = [mdl_vec(1:end-2), '];'];
	mdl_vec = eval(mdl_vec);
	
	[p, rnorm, r, exitflag, output, lambda, J] = lsqcurvefit(mdl_vec, p0, x_vec, y_vec, lb, ub, varargin{1:end});
end