function [ax] = barplot(ax, xpos, yvals, options, barProperties, axProperties)
	% Generates a bar chart w/ error bars and (optionally) individual points)
	%
	%   Inputs:
	%
	%		ax				<axes> Axes handle to plot on
	%
	%		xpos			<numeric> A vector of x-axis centers for each bar
	%
	%       yvals			<numeric> A matrix of values for each point
	%							Note: mean and standard deviation are
	%							automatically computed to plot the bar and
	%							errorbars. To change to geomean or geostd,
	%							use the 'options' input. 
	%
	%		options			<cell, char> (Optional) A cell array of strings (or 
	%						a single string) specifying optional plotting behavior:
	%							'geomean'	Compute bar height w/ geomean 
	%										rather than mean.
	%							'geostd'	Compute errorbar height w/
	%										geostd rather than std.
	%							'points'	Show individual points
	%							'unierr'	Plot errorbars only on the
	%										side of the bar facing away
	%										from the x-axis
	%
	%       barProperties   <cell> (Optional) Cell array of bar plot optional inputs
	%
	%		axProperties	<struct> (Optional) Struct of axes properties
	%
	%   Outputs: 
	%
	%       ax              A handle to the figure axes
	% 
	% Written By 
	% Ross Jones
	% jonesr18@mit.edu
	% Weiss Lab, MIT
	%
	% Update Log:
	%   

	[ymeans, ystdsNeg, ystdsPos, errProperties] = zCheckInputs_barPlot();

	hold(ax, 'on')
	bar(ax, xpos, ymeans, barProperties{:})
	errorbar(ax, xpos, ymeans, ystdsNeg, ystdsPos, '.', 'markersize', 1, errProperties{:})

	if any(ismember({'points', 'dots'}, options))
		% Copmute an average dx so that we know how much space there is
		% between each bar for plotting dots
		dx = mean(diff(xpos));

		xvals = linspace(0.7, 1.3, size(yvals, 1))';
		rxi = randperm(numel(xvals))';

		xvals = repmat(xpos, size(yvals, 1), 1) + dx .* (xvals(rxi) - 1);

		scatter(ax, xvals(:), yvals(:), 10, 'k', 'filled')
	end

	% Set axes properties
	for f = fieldnames(axProperties)'
		if isstruct(axProperties.(f{:}))
			set(ax.(f{:}), axProperties.(f{:}))
		else
			set(ax, f{:}, axProperties.(f{:}));
		end
	end


	% --- Helper Functions --- %


	function [ymeans, ystdsNeg, ystdsPos, errProperties] = zCheckInputs_barPlot()

		if ~exist('barProperties', 'var')
			barProperties = {};
		end
		if ~exist('axProperties', 'var') || isempty(axProperties)
			axProperties = struct();
		end
		if exist('options', 'var')
			if ischar(options), options = {options}; end
		else
			options = {};
		end

		if ismember('geomean', options)
			ymeans = 10.^mean(log10(yvals), 1, 'omitnan');
		else
			ymeans = mean(yvals, 1, 'omitnan');
		end

		if ismember('geostd', options)
			ystds = geostd(yvals, 0, 1, 'omitnan');
		else
			ystds = std(yvals, 0, 1, 'omitnan');
		end

		if ismember ('unierr', options)
			% Make errors only show up on the appropriate 'side' of the bar
			ystdsNeg = ystds .* (ymeans < 0);
			ystdsPos = ystds .* (ymeans > 0);
		else
			ystdsNeg = ystds;
			ystdsPos = ystds;
		end


		eci = find(strcmpi('edgecolor', barProperties));
		if ~isempty(eci)
			errProperties = {'color', barProperties{eci + 1}};
		else
			errProperties = {'color', 'k'};
			barProperties = [barProperties, {'edgecolor', 'k'}];
		end

		lci = find(strcmpi('linewidth', barProperties));
		if ~isempty(lci)
			errProperties = [errProperties, {'linewidth', barProperties{lci + 1}}];
		else
			errProperties = [errProperties, {'linewidth', 1}];
			barProperties = [barProperties, {'linewidth', 1}];
		end

	end

end