function [relErrNeg, relErrPos] = relError(means, stdevs)
	% Computes 'relative error': used when plotting errorbars on log plots
	%
	%	[relErrNeg, relErrPos] = relError(means, stdevs)
	%
	%	Inputs
	%
	%		means				<numeric> A matrix of mean values
	%
	%		stdevs				<numeric> A matrix of standard deviation values
	%
	%	Outputs
	%
	%		relErrNeg			<numeric> Negative errorbar magnitude
	%
	%		relErrPos			<numeric> Positive errorbar magnitude
	%
	%	Example
	%
	%		means = mean(ydata, 2, 'omitnan');
	%		stdevs = std(ydata, 2, 'omitnan');
	%		[relErrNeg, relErrPos] = relError(means, stdevs);
	%		figure(), ax = gca();
	%		errorbar(ax, xdata, means, relErrNeg, relErrPos)
	%		set(ax, 'YScale', 'log')
	%
	% Written By
	% Ross Jones
	% jonesr18@mit.edu
	% Weiss Lab, MIT
	
	relErr = 1 / log(10) * stdevs ./ means;
	relErrPos = 10.^(log10(means) + relErr) - means;
	relErrNeg = means - 10.^(log10(means) - relErr);

end