classdef Plotting < handle
	% Compilation of plotting functions from the Weiss Lab flow cytometry MATLAB repository
	%
	%	Methods are implemented as static functions, so they can be called directly
	%	or after creating a handle to this class.
	%
	%		Ex:
    %       ax = gca();
	%		Plotting.biexpAxes(ax);
	%
	%		Ex2:
	%		plt = Plotting();
	%		flow.densityplot(ax, xdata, ydata, nPoints, mode, colorMap)
    %
    %   Functions:
    %
    %       densityplot(ax, xdata, ydata, nPoints, mode, colorMap)
    %       violinplot(ax, ydata, xcenter, mode, faceColor)
    %       biexplot(x, y, plotArgs, options)
    %       biexpAxes(ax, biexpX, biexpY)
    %       [nelements , centers] = biexhist(Y, M)
    %       handle = histFit(data, channel, dataType, faceColor, axPosition)
    %       scatterBarPlot(ax, data, channel, dataType, colorMaps, transformation)
    %       batchViolinPlot(ax, data, channel, dataType, colors)
    %       lineDensityPlot(struct('ax', 'data', 'channel', 'dataType', 'colors', 'xscale', 'shade'))
    %       orthogonalityMatrix(A,labels)
    %       [ax, h] = standardHeatmap(data, cmap, rowLabels, colLabels, norm, doCluster, symmetric)
    %
    % Written/Compiled by
	% Ross Jones
	% jonesr18@mit.edu
    % Weiss Lab, MIT
    
	methods (Static)
		
		
		function ax = setAxProps(ax, axProperties)
			% An easier method to apply a struct of axes properties to an axis
			%
			%	ax = Plotting.setAxProps(ax, axProperties);
			%
			%	Inputs
			%
			%		ax				<Axes> The axes object to modify
			%
			%		axProperties	<struct> A struct array of axes property-value
			%						pairs to modify. 
			%
			%	Outputs
			%
			%		ax				<Axes> The updated axes object
			%
			%	Axes handles point to objects with many fields, some of which
			%	are themselves struct arrays. To edit these structs, we cannot
			%	just say set(ax, 'Title', struct( ... )) - the Axes gets
			%	confused because it is expecting a full String with several
			%	other fields as well. Instead, we have to edit each desired
			%	property independently. This method handles that while allowing
			%	the input to still be a simple struct of property-value pairs. 
			
			for f = fieldnames(axProperties)'
				if isstruct(axProperties.(f{:}))
					set(ax.(f{:}), axProperties.(f{:}))
				else
					set(ax, f{:}, axProperties.(f{:}));
				end
			end
		end
		
		
		function [density, binCenters, binCounts] = computeDensity(dataMatrix, dMode, numPoints, nonZero)
			% Computes the density of a given set of points
			%
			%	density = computeDensity(dataMatrix, mode, numPoints, nonZero)
			%
			%	Inputs
			%		dataMatrix		<numeric> An NxM matrix of N elements in M 
			%						dimensions for which to compute density
            %
            %       dMode			<char> Determines how density is calculated
            %                         'neighbors'     Based on number of nearby points
            %                         'fastn'         'neighbors' using 1/10 points for speed
            %                         'kernel'        Kernel density estimation 
            %                                         (auto-selects bandwidth - see kde2d.m)
            %                         'hist'          Interpolates density from a 2D histogram
            %       
			%		numPoints		(Optional) <numeric> The number of elements to use 
			%						for the density calculation (for speedup)
			%
			%		nonZero			(Optional) <logical> Flag to only estimate using
			%						non-zero values
			%
			%	Outputs
			%
			%		density			An Nx1 vector of density values for each
			%						element in dataMatrix.
			%
			%		binCenters		If computing with 'hist' mode, then the bin centers 
			%						used in the calculation can be requested. 
			%
			%		binCounts		As with above, the binCounts can also be requested
			%
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
			
			% Check inputs
			dataValid = zCheckInputs_computeDensity();
            numDims = size(dataMatrix, 2);
			
			% Check mode of point coloration
			subsample = randperm(size(dataValid, 1), numPoints);
			dataSub = dataValid(subsample, :);
			binCenters = [];	% Only fill if using dMode = 'hist'
			binCounts = [];		% Same as above
			switch dMode 
				case {'neighbors', 'fastn'}
					% Find how many neighbours there are less than dX and dY away.
					if (strcmpi(dMode, 'fastn'))
						skip = 10;
					else
						skip = 1;
					end
					
					% Find the range of each dimension
					% Use valid data for this mode because we aren't going to
					% return the subsample vector
					dist = max(dataValid, 1) - min(dataValid, 1);
					
					% Find the density by dividing the number of points by the
					% range. Dividing 100 by this number defines some distance  
					% one would expect close data points to be from one another.
					dens = numPoints ./ dist;
					dX = 100 ./ dens;
					
					neighbors = zeros(numPoints, 1);
					for p = 1:size(dataValid, 1) % Don't use numPoints
						vals = dataValid(p, :);
						
						closePoints = true(ceil(numPoints / skip), 1);
						for d = 1:numDims
							closePoints = closePoints & ...
									(abs(dataSub(1:skip:end, d) - vals(d)) < dX(d));
						end
						neighbors(p) = sum(closePoints) * skip;
					end
					
%					% Convert density to log scale to get better view of data
% 		            interpDensity = log10(interpDensity);
					
					density = neighbors;

				case 'kernel'
					
					if (numDims > 2)
						warning('Kernel Density only valid for <= 2 dimensions')
					end
					
					% Estimate density with kernel
					if (numDims == 1)
						[bw, kdensity, meshX] = kde(dataSub);
						interpDensity = interp1(meshX, kdensity, dataMatrix);
					else
						[bw, kdensity, meshX, meshY] = kde2d(dataSub(:, [1, 2]));
						interpDensity = interp2(meshX, meshY, kdensity, dataMatrix(:, 1), dataMatrix(:, 2));
					end
					fprintf(1, 'Kernal density estimated with bandwidth: %.3f\n', bw);
					
% 		            % Convert density to log scale to get better view of data
% 		            interpDensity = log10(interpDensity);
					
					density = interpDensity;
					
				case 'hist'
					
					numBins = 25;
					binEdges = cell(1, numDims);
					for d = 1:numDims
						% Add a little extra padding to ensure all data is
						% captured. Was having issues with small values getting
						% chopped off and density to be returned as NaN
						sh = range(dataValid(:, d)) / 20;
						binEdges{d} = linspace(min(dataValid(:, d)) - sh, max(dataValid(:, d)) + sh, numBins + 1);
					end
					
					% This is actually the best way to do N-dimensional
					% histogram computation - it's just like binning!
					dataInBins = FlowAnalysis.simpleBin(dataSub, binEdges);
					binCounts = cellfun(@numel, dataInBins);
					
					% Estimate bin centers by averaging the edges
					binCenters = cell(1, numDims);
					for d = 1:numDims
						binCenters{d} = zeros(1, numBins);
						for b = 1:numBins
							binCenters{d}(b) = mean(binEdges{d}([b, b + 1]));
						end
					end
					
					% Create grid (1xnumDims cell array)
					grid = createGrid(binCenters);
					
					% Create cell version of data so we can throw it into
					% interpn without knowing how many dimensions there are
					dataCell = cell(1, numDims);
					for d = 1:numDims, dataCell{d} = dataMatrix(:, d);	end
					
% 					sizes = cellfun(@size, grid, 'uniformoutput', false);
% 					sizes{:}
% 					size(binCounts)
					
					% Interpolate over the 2D histogram counts to make a PDF for the data points
					interpDensity = interpn(grid{:}, binCounts, dataCell{:});
					
% 		            % Convert density to log scale to get better view of data
% 		            interpDensity = log10(interpDensity);
					
					density = interpDensity;
					
			end
			
			
			% --- Helper Functions --- %
			
			
			function dataValid = zCheckInputs_computeDensity()
				
				validateattributes(dataMatrix, {'numeric'}, {}, mfilename, 'data', 1);
				validatestring(dMode, {'neighbors', 'fastn', 'kernel', 'hist'}, mfilename, 'mode', 2);
				
				if exist('nonZero', 'var')
					nonZero = any(logical(nonZero(:)));
					validateattributes(nonZero, {'logical'}, {}, mfilename, 'nonZero', 4);
				else
					nonZero = false;
				end
				
				 % Get absolute value of complex values
				dataMatrix = sign(dataMatrix) .* abs(dataMatrix);

				% Ignore NaN and inf values
				valid = ~(sum(isnan(dataMatrix), 2) | sum(isinf(dataMatrix), 2));

				% Ignore negative values if requested
				if (nonZero), valid = (valid & (sum(dataMatrix >= 0, 2) > numDims)); end

				% Only need to check x and y for valid points
				dataValid = dataMatrix(valid, :);
				
				if exist('numPoints', 'var')
					validateattributes(numPoints, {'numeric'}, {'scalar', 'positive'}, mfilename, 'numPoints', 3);
				else
					numPoints = inf;
				end
				numPoints = round(min(size(dataValid, 1), numPoints));

			end
			
			
			function grid = createGrid(edges)
				% Create a grid marking the binned space
				% --> Based on ndgrid.m, which I can't use becuase it
				%	  works based on variable outputs.
				
				gridSize = cellfun(@numel, edges);
				
				grid = cell(1, numel(edges));
				for g = 1:numel(edges)
					% Reshape edges into the desired dimension
					x = edges{g};
					s = ones(1, numel(edges{g}) - 1);
					s(g) = numel(x);
					x = reshape(x, s);

					% Repeat edges into every other dimension except d
					gs = gridSize; 
					gs(g) = 1;
					grid{g} = repmat(x, gs);
				end
			end
		end
		
		
		function [colors, sortIdx] = getColors(inputData, cmap, options)
			% Generates a color representing the value of each element of a
			% given input vector/matrix using a given color scheme
			%
			%	[colors, sortIdx] = getColors(inputData, cmap, options)
			%
			%	Inputs
			%		inputData	<numeric> A vector/matrix of data values for 
			%					which to generate colors.
			%
			%		cmap		<numeric, char, ColorMap> (Optional) 
			%					The colormap to represent data values. 
			%					 - Can input an Nx3 matrix of RGB values, a
			%					   ColorMap object pre-initialized with a
			%					   color, or a string indicating which
			%					   ColorMap to initialize. 
			%					 - ColorMap and char inputs will yield 100
			%					   unique color values on the given scale
			%					 - Defualt = parula(100);
			%
			%		options		<struct> A number of optional inputs:
			%						'min': The minimum color threshold 
			%							(default = min(inputData))
			%						'max': The maximum color threshold 
			%							(default = max(inputData))
			%						'nan': The color to set NaN values to. 
			%						Accepted options: 'white', 'grey', 'black', <color>
			%						where <color> is an RGB triplet
			%							(default = 'white')
			%						'imag': The color to set complex values to. 
			%						Accepted options: 'white', 'grey', 'black', <color>
			%						where <color> is an RGB triplet
			%							(default = 'grey')
			%
			%	Outputs
			%
			%		colors		An NxMx...x3 matrix of color values corresponding
			%					with the input data values, where the appended
			%					dimension is in the order RGB. Row-vector inputs
			%					will return RGB values in independent rows. 
			%					*** Returned unsorted unless 'sortIdx' is
			%						requested as an output!
			%					*** NaN and complex inputData values are
			%						converted to white. If this behavior is not
			%						desired, change the 'nan' and 'imag' options
			%
			%		sortIdx		The sorting order for inputData. 
			%					*** If this output is requested, then 'colors'
			%						is returned in sorted order!
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
			
			% Check inputs
			[inputSize, MIN, MAX, nanInfo, imagInfo] = zCheckInputs_getColors();
			numColors = size(cmap, 1);
			
			% Convert values to color
			colorIdxs = max(min(ceil((inputData - MIN) ./ (MAX - MIN) .* numColors), numColors), 1);
			try
				colors = cmap(colorIdxs, :);
			catch ME
				fprintf(2, 'Min: %d\n', MIN)
				fprintf(2, 'Max: %d\n', MAX)
				fprintf(2, 'Min index: %d\n', min(colorIdxs))
				fprintf(2, 'Max index: %d\n', max(colorIdxs))
				fprintf(2, 'NaNs present: %d\n', any(isnan(colorIdxs)))
				fprintf(2, 'All vals real: %d\n', all(isreal(colorIdxs)))
				error('Encountered error getting colors')
				
			end
			
			% Set NaN/Imag colors
			nanColors = repmat(nanInfo.color, size(colors, 1), 1);
			colors(nanInfo.idxs, :) = nanColors(nanInfo.idxs, :);
			imagColors = repmat(imagInfo.color, size(colors, 1), 1);
			colors(imagInfo.idxs, :) = imagColors(imagInfo.idxs, :);
			
			% If sortIdx output requested, sort so brightest values come out last
			if nargout > 1
				[~, sortIdx] = sort(inputData);
				colors = colors(sortIdx, :);
			end
			
			% Re-shape back to the original input dimensions, and put 
			% the colors in the N+1th dimension.
			if (inputSize(2) == 1)
				reshapeSize = [inputSize(1), 3];	% Handles column vector input
			elseif (inputSize(1) == 1)
				reshapeSize = [inputSize(2), 3];	% Handles row vector input
			else
				reshapeSize = [inputSize, 3];		% Handles 2D+ input
			end
			colors = reshape(colors, reshapeSize);
			
			
			% --- Helper Functions --- %
			
			
			function [inputSize, MIN, MAX, nanInfo, imagInfo] = zCheckInputs_getColors()
				
				validateattributes(inputData, {'numeric'}, {}, mfilename, 'inputData', 1);
				inputSize = size(inputData);
				inputData = reshape(inputData, [], 1); % For simplicity
				
				if exist('cmap', 'var')
					cmap = Plotting.checkCmap(cmap);
				else
					cmap = parula(100);
				end
								
				% Treat NaNs (no data) and complex (log(neg values)) as min 
				% values, then later overwrite them with white so they don't 
				% show up on normal plots
				nanIdxs = isnan(inputData);
				inputData(nanIdxs) = min(inputData(~isinf(inputData)));
				imagIdxs = (imag(inputData) ~= 0);
				inputData(imagIdxs) = min(inputData(~imagIdxs & ~isinf(inputData)));
				
				% Default options
				MIN = min(inputData(~isinf(inputData)));
				MAX = max(inputData(~isinf(inputData)));
				nanColor = [1, 1, 1];	
				imagColor = [0.5, 0.5, 0.5];
				
				if exist('options', 'var')
					if isfield(options, 'min')
						MIN = options.min;
						validateattributes(MIN, {'numeric'}, {'scalar'}, mfilename, 'options.min', 3)
					end
					
					if isfield(options, 'max')
						MAX = options.max;
						validateattributes(MAX, {'numeric'}, {'scalar'}, mfilename, 'options.max', 3)
					end
					
					if isfield(options, 'nan')
						nanColor = options.nan;
						validateattributes(nanColor, {'numeric', 'char'}, {}, mfilename, 'options.nan', 3)
						if ischar(nanColor)
							if strcmpi(nanColor, 'white')
								nanColor = [1, 1, 1];
							elseif strcmpi(nanColor, 'grey')
								nanColor = [0.5, 0.5, 0.5];
							elseif strcmpi(nanColor, 'black')
								nanColor = [0, 0, 0];
							else
								error('NaN color not recognized')
							end
						else
							nanColor = reshape(nanColor, 1, []);
							assert(numel(nanColor) == 3);
							if (min(nanColor) < 0), nanColor = nanColor + min(nanColor); end
							if (max(nanColor) > 1), nanColor = nanColor ./ max(nanColor); end
						end
					end
					
					if isfield(options, 'imag')
						imagColor = options.imag;
						validateattributes(imagColor, {'numeric', 'char'}, {}, mfilename, 'options.imag', 3)
						if ischar(imagColor)
							if strcmpi(imagColor, 'white')
								imagColor = [1, 1, 1];
							elseif strcmpi(imagColor, 'grey')
								imagColor = [0.5, 0.5, 0.5];
							elseif strcmpi(imagColor, 'black')
								imagColor = [0, 0, 0];
							else
								error('Imag color not recognized')
							end
						else
							imagColor = reshape(imagColor, 1, []);
							assert(numel(imagColor) == 3);
							if (min(imagColor) < 0), imagColor = imagColor + min(imagColor); end
							if (max(imagColor) > 1), imagColor = imagColor ./ max(imagColor); end
						end
					end
				end
				
				% Fix for inf values
				inputData(inputData == inf) = MAX;
				inputData(inputData == -inf) = MIN;
				
				nanInfo = struct('idxs', nanIdxs, 'color', nanColor);
				imagInfo = struct('idxs', imagIdxs, 'color', imagColor);
			end
		end
		
		
		function densityplot(ax, xdata, ydata, nPoints, dMode, cmap, dotSize, nonZero)
			% Plots a scatterplot with colors of the dots indicating density
			%
			%	densityplot(ax, xdata, ydata, nPoints, mode, cmap, dotSize, nonZero)
			%   
            %   Inputs
            %   
			%       xdata (vector)      X-axis values
			%       
            %       ydata (vector)      Y-axis values the same size as xdata
            %
            %       nPoints (integer)   The number of points to plot
            %                            - automatically selects min(nPoints, length(xdata))
            %
            %       dMode (string)		Determines how density is calculated
            %							 'neighbors'     Based on number of nearby points
            %							 'fastn'         'neighbors' using 1/10 points for speed
            %							 'kernel'        Kernel density estimation 
            %								             (auto-selects bandwidth - see kde2d.m)
            %							 'hist'          Interpolates density from a 2D histogram
            %                            <numerical>     Colors the points based on a given set of values
            %                                            Must be the same size as xdata/ydata
            %       
			%		cmap (numeric, char, ColorMap) (Optional) 
			%							The colormap to represent data values. 
			%							 - Can input an Nx3 matrix of RGB values, a
			%							   ColorMap object pre-initialized with a
			%							   color, or a string indicating which
			%							   ColorMap to initialize. 
			%							 - ColorMap and char inputs will yield 100
			%							   unique color values on the given scale
			%							 - Defualt = parula(100);
			%
			%		dotSize (integer)	(Optional) The size of dots to use (default = 8)
			%
			%		nonZero (logicle)	(Optional) Indicates whether to force
			%							non-zero values
            %
			%
			%   Example:
			%       [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs('sample.fcs');
			%       GreenChannel = getChannel(fcshdr,'FIT');
			%       RedChannel = getChannel(fcshdr,'Red');
			%       GreenData = fcsdat(:,GreenChannel);
			%       RedData = fcsdat(:,RedChannel);
            %       figure();
            %       ax = subplot(1, 5, 1);
			%       densityplot(ax, greenData, redData, 5000, 'kernel', ColorMap('red'));
			%
			% Written By
			% Breanna DiAndreth
			% bstillo@mit.edu
			% Weiss Lab, MIT
			%
			% Update Log:
			%	2015-02-06 (Ross) - Doubled speed by slightly changing how neighbours is calculated
            %	2016-03-28 (Ross) - Added kernel density estimation and merged with fastdensity
            %                          calculation function
            
            % Check inputs
            zCheckInputs_densityplot();
			
			% Determine if density-based or directly-supplied coloration
			if ischar(dMode)
				density = Plotting.computeDensity([xdata, ydata], dMode, nPoints, nonZero);
				[colors, sortIdx] = Plotting.getColors(density, cmap);
			else
				[colors, sortIdx] = Plotting.getColors(dMode, cmap, struct('min', 0, 'max', 4.5, 'numColors', 100));
			end
			
			ss = FlowAnalysis.subSample(numel(xdata), nPoints);
			scatter(ax, xdata(sortIdx(ss)), ydata(sortIdx(ss)), dotSize, colors(ss, :), 'filled');
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_densityplot()
				assert(all(size(xdata) == size(ydata)), 'xdata and ydata are different sizes!')
				validateattributes(xdata, {'numeric'}, {'vector'}, mfilename, 'xdata', 2);
				validateattributes(ydata, {'numeric'}, {'vector'}, mfilename, 'ydata', 3);
				validateattributes(nPoints, {'numeric'}, {'scalar'}, mfilename, 'nPoints', 4);
				nPoints = round(nPoints);
				validateattributes(dMode, {'char', 'numeric'}, {'vector'}, mfilename, 'mode', 5);
				if ~ischar(dMode)
					assert(all(size(dMode) == size(xdata)), 'If numeric, mode must be the same size as the data!')
				end
				
				% Check colormap (bulk of checking happends in getColors())
				if exist('cmap', 'var')
					validateattributes(cmap, {'numeric', 'char', 'ColorMap'}, {}, mfilename, 'cmap', 6);
				else
					cmap = parula(100);
				end
				
				if exist('dotSize', 'var')
					validateattributes(dotSize, {'numeric'}, {'scalar'}, mfilename, 'dotSize', 7);
% 					dotSize = round(dotSize);
				else
					dotSize = 8;
				end
				
				if exist('nonZero', 'var')
					nonZero = any(logical(nonZero(:)));
					validateattributes(nonZero, {'logical'}, {}, mfilename, 'nonZero', 8);
				else
					nonZero = false;
				end
			end
		end
        
        
		function h = violinplot(ax, ydata, xpos, binEdges, faceColor)
			% Plots the data as a violin plot, which represents density as a vertical,
			% horizontally symmetric histogram. The histograms are overlaid with
			% a black box-and-whisker plot with the box representing the 25th - 75th
			% percentiles, the whiskers showing the 95th - 5th percentiles, and the
			% white dot in the middle showing the median.
			%   
			%	violinplot(ax, ydata, xpos, binEdges, faceColor)
			%
			%   Inputs
			%   
			%		ax (handle)			The axis handle to plot to
			%
			%       ydata (vector)      Values of data to plot
			%
			%       xpos (integer)		The x-position of the violin
			%
			%		binEdges (vector)	Bin edges for histogram computation
			%							 - Should be equal for all violins for
			%							   valid comparison
			%       
			%       faceColor (rgb)     The color of the violin face. The face is plotted with an
			%                           alpha value of 0.5 to look nice.
			%
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			%
			% Update Log:
			%
			%	2018-01-28:		Switched to using separate density computing function
			%	2018-10-03:		Fixed density calculation and added box plot overlay

			% Check inputs
			zCheckInputs_violinplot(); 
			
			% Compute density
			binCounts = histcounts(ydata, binEdges, 'Normalization', 'count');
			binCenters = reshape(binEdges(2:end) - diff(binEdges) / 2, [], 1);
			
			% Normalize to max value, then multiply by 0.4 so that 
			%	violins don't overlap
% 			binCounts = reshape(log10(binCounts), [], 1) * 0.4 / max(log10(binCounts(:)));
			binCounts = reshape((binCounts), [], 1) * 0.4 / max((binCounts(:)));
			valid = ~isinf(binCounts) & ~isnan(binCounts);
			
			% Plot violins
			h = fill(ax, [xpos + binCounts(valid); xpos - flipud(binCounts(valid))], ...
					[binCenters(valid); flipud(binCenters(valid))], ...
					faceColor, 'EdgeColor', 'none', 'LineWidth', 1);
			
			% Plot pseudo boxplot
			prctiles = prctile(ydata, [95, 75, 50, 25, 5]);
			errorbar(ax, xpos, prctiles(3), prctiles(3) - prctiles(5), ...
					prctiles(1) - prctiles(3), 'k', 'linewidth', 0.5)
			fill(ax, [xpos + 0.07 * ones(2, 1); xpos - 0.07 * ones(2, 1)], ...
					[prctiles([2, 4])'; prctiles([4, 2])'], ...
					'k', 'EdgeColor', 'none');
			plot(ax, xpos, prctiles(3), '.w', 'markersize', 16);
			
			% TODO Get this to work someday - technically better than the hacked
			% version above since the whiskers should be calculated by an algorithm
% 			boxplot(ax, ydata, 'positions', xcenter, 'boxstyle', 'filled', 'colors', 'k')
			
			
			% --- Helper Function --- %


			function zCheckInputs_violinplot()
				
				validateattributes(ax, {'matlab.graphics.axis.Axes'}, {}, mfilename, 'ax', 1);
				
				validateattributes(ydata, {'numeric'}, {'vector'}, mfilename, 'ydata', 2);
				ydata = reshape(ydata, [], 1); % Force column vector
				
				% Fix negative infinite values by setting the resulting values to the 
				% otherwise minimum value.
				if any(ydata == -inf)
					warning('Negative values detected - setting to min value')
					ydata(ydata == -inf) = min(ydata(ydata ~= -inf));
				end
				
				validateattributes(xpos, {'numeric'}, {}, mfilename, 'xpos', 3);
				xpos = round(xpos(1)); % Ensure integer and only one point
				
				validateattributes(binEdges, {'numeric'}, {'vector'}, mfilename, 'binEdges', 4);
				validateattributes(faceColor, {'numeric'}, {'vector'}, mfilename, 'faceColor', 5);
				assert(length(faceColor) == 3);
				
			end
		end
        
		
		function h = violincompare(ax, ydata1, ydata2, xcenter, binEdges, faceColors, options)
			% Plots the data as a split violin plot comparing two samples.
			% The histograms are overlaid with lines representing the 5th, 
			% 25th, 50th (median), 75th, and 95th percentiles.
			%   
			%	violincompare(ax, ydata, xcenter, binEdges, faceColor)
			%
			%   Inputs
			%   
			%		ax (handle)			The axis handle to plot to
			%
			%       ydata1 (vector)		Values of data to plot 
			%							(Sample 1, plotted on the left)
			%
			%       ydata2 (vector)		Values of data to plot
			%							(Sample 2, plotted on the left)
			%
			%       xcenter (integer)   The value to center the x-values on
			%
			%		binEdges (vector)	Bin edges for histogram computation
			%							 - Should be equal for all violins for
			%							   valid comparison
			%       
			%       faceColors (rgb)	2x3 matrix of colors for violin faces.
			%							Row 1: left face (ydata1)
			%							Row 2: right face (ydata2)
			%							The faces are plotted with an alpha value 
			%							of 0.5 to look nice.
			%
			%		options (cell/char) (optional) A cell array of strings (or a single 
			%							string) specifying optional plotting behavior:
			%							'smooth'	Flag to smooth the histogram data
			%
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			%
			% Update Log:
			%
			%	2018-01-28:		Switched to using separate density computing function
			%	2018-10-03:		Fixed density calculation and added box plot overlay

			% Check inputs
			zCheckInputs_violincompare();
			
			% Compute density
			binCounts1 = histcounts(ydata1, binEdges, 'Normalization', 'count');
			binCounts2 = histcounts(ydata2, binEdges, 'Normalization', 'count');
			binCenters = binEdges(2:end) - diff(binEdges) / 2;
			
			% Normalize to max value, then multiply by 0.4 so that 
			%	violins don't overlap
% 			binCounts = log10(binCounts) * 0.4 / max(log10(binCounts));
			if ismember('smooth', options)
				binCounts1 = smooth(binCounts1);
				binCounts2 = smooth(binCounts2);
			end
			binCounts1 = binCounts1 * 0.4 / max(binCounts1(:));
			binCounts2 = binCounts2 * 0.4 / max(binCounts2(:));
			
			% Force first and last values to be at zero to avoid messed-up plots
			binCounts1 = [0; binCounts1; 0];
			binCounts2 = [0; binCounts2; 0];
			binCenters = binCenters([1, 1:end, end]);
			
			% Throw out inf and nan values
			valid1 = ~isinf(binCounts1) & ~isnan(binCounts1);
			valid2 = ~isinf(binCounts2) & ~isnan(binCounts2);

			% Plot violins (ydata1 on left, ydata2 on right)
			h = zeros(1, 2);
			h(1) = fill(ax, xcenter - binCounts1(valid1), binCenters(valid1), ...
					faceColors(1, :), 'EdgeColor', 'none', 'LineWidth', 1);
			h(2) = fill(ax, xcenter + binCounts2(valid2), binCenters(valid2), ...
					faceColors(2, :), 'EdgeColor', 'none', 'LineWidth', 1);
			
			% Plot pseudo boxplot
			pcts1 = prctile(ydata1, [95, 75, 50, 25, 5]);
			pcts2 = prctile(ydata2, [95, 75, 50, 25, 5]);
			
			linewidths = [0.5, 1, 2, 1, 0.5];
			colorSums = mean(faceColors, 2);
% 			lineColors = double((colorSums < 0.4) .* ones(2, 3));
			lineColors = double((colorSums < 0.4) .* zeros(2, 3));
			for pi = 1:numel(linewidths)
				[~, ci1] = min(abs(binCenters - pcts1(pi)));
				[~, ci2] = min(abs(binCenters - pcts2(pi)));
				
				plot(ax, xcenter - [binCounts1(ci1), 0], binCenters([ci1, ci1]), ...
						'color', lineColors(1, :), 'linewidth', linewidths(pi))
				plot(ax, xcenter + [0, binCounts2(ci2)], binCenters([ci2, ci2]), ...
						'color', lineColors(2, :), 'linewidth', linewidths(pi))
			end
			
			
			% --- Helper Function --- %


			function zCheckInputs_violincompare()

				validateattributes(ydata1, {'numeric'}, {'vector'}, mfilename, 'ydata1', 2);
				validateattributes(ydata2, {'numeric'}, {'vector'}, mfilename, 'ydata2', 3);
				ydata1 = reshape(ydata1, [], 1); % Force column vector
				ydata2 = reshape(ydata2, [], 1); % Force column vector
				
				% Fix negative infinite values by setting the resulting values to the 
				% otherwise minimum value.
				if any(ydata1 == -inf)
					warning('Negative values detected - setting to min value')
					ydata1(ydata1 == -inf) = min(ydata1(ydata1 ~= -inf));
				end
				if any(ydata2 == -inf)
					warning('Negative values detected - setting to min value')
					ydata2(ydata2 == -inf) = min(ydata2(ydata2 ~= -inf));
				end
				
				validateattributes(xcenter, {'numeric'}, {}, mfilename, 'xcenter', 4);
				if (round(xcenter) ~= xcenter)
					warning('Non-integer ''xcetner'', rounding');
					xcenter = round(xcenter(1)); % Ensure integer and only one point
				end
				
				validateattributes(binEdges, {'numeric'}, {'vector'}, mfilename, 'binEdges', 5);
				binEdges = reshape(binEdges, [], 1); % Force column vector
				validateattributes(faceColors, {'numeric'}, {}, mfilename, 'faceColor', 6);
				assert(all(size(faceColors) == [2, 3]), ...
						'Input ''faceColors'' must be a 2 x 3 matrix');
				
				if ~exist('options', 'var')
					options = {};
				end
			end
		end
		
        
		function biexplot(x, y, plotArgs, options)
			% NOTE: #Decrepit - to be removed before full release
			%
			%BIEXPLOT(...) is the same as PLOT(...) except that logicle scales (Parks,
			%et al.) are used for both the X- and Y- axes.
			%   As with PLOT(...), various line types, plot symbols and colors may be obtained with
			%   BIEXPLOT(X,Y,S) where S is a character string made from one element
			%   from any or all the following 3 columns:
			%
			%          b     blue          .     point              -     solid
			%          g     green         o     circle             :     dotted
			%          r     red           x     x-mark             -.    dashdot 
			%          c     cyan          +     plus               --    dashed   
			%          m     magenta       *     star             (none)  no line
			%          y     yellow        s     square
			%          k     black         d     diamond
			%          w     white         v     triangle (down)
			%                              ^     triangle (up)
			%                              <     triangle (left)
			%                              >     triangle (right)
			%                              p     pentagram
			%                              h     hexagram
			%
			%   When the string 'density' is supplied for S, a density plot is created
			%   using logicle scales.  This most closely imitates FlowJo plots
			%
			%   Example:
			%       [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs('sample.fcs');
			%       GreenChannel = getChannel(fcshdr,'FIT');
			%       RedChannel = getChannel(fcshdr,'Red');
			%       GreenData = fcsdat(:,GreenChannel);
			%       RedData = fcsdat(:,RedChannel);
			%       biexplot(greenData,redData,'density')
			%
			% Written By
			% Breanna DiAndreth
			% bstillo@mit.edu
			%
			% Update Log:
			%	2016-02-22 (Ross):		Cleaned up some code
            
			% Check inputs
			validateattributes(x, {'numeric'}, {}, mfilename, 'x', 1);
			validateattributes(y, {'numeric'}, {}, mfilename, 'y', 2);
			if (exist('plotArgs', 'var'))
				validateattributes(plotArgs, {'cell', 'char'}, {}, mfilename, 'plotArgs', 3);
			end
			if (exist('options', 'var'))
				validateattributes(options, {'numeric'}, {}, mfilename, 'options', 1);
			end
			
			xlog = Transforms.lin2logicle(x);
			ylog = Transforms.lin2logicle(y);
			
			if (~exist('plotArgs', 'var'))
				plotArgs = {'.', 'MarkerSize', 2};
				plot(xlog, ylog, plotArgs{:})
			else
				if (any(strcmpi(plotArgs, 'density')))
					Plotting.densityplot(gca(), xlog, ylog, numel(xlog), 'normal', ColorMap('parula'))
                elseif (any(strcmpi(plotArgs, 'fastdensity')))
                    Plotting.densityplot(gca(), xlog, ylog, 5000, 'fast', ColorMap('parula'))
                elseif (any(strcmpi(plotArgs, 'kerneldensity')))
                    Plotting.densityplot(gca(), xlog, ylog, 5000, 'kernel', ColorMap('parula'))
                else
					plot(xlog, ylog, plotArgs{:})
				end
			end
			
			if (exist('options', 'var'))
				if (options == 1)
		%             fprintf('Axis correction... may take a little time...please wait.\n')
					% --> not anymore!
					hha1 = gca;
                    Plotting.biexpAxes(hha1);
				end
			end
        end
		
        
        function ax = biexpAxes(ax, biexpX, biexpY, biexpZ, doMEF, params)
            % Tranforms the given axes to a scale used for biexponential views.
            %   
			%	ax = biexpAxes(ax, biexpX, biexpY, biexpZ, doMEF, params);
			%
            %   Inputs (all optional - defaults to X/Y both logicle)
            %       
            %       biexpX	(logical)		TRUE - make x-axis biexponential
            %       biexpY	(logical)       TRUE - make y-axis biexponential
            %       biexpZ	(logical)       TRUE - make z-axis biexponential
			%		doMEF	(logical)		TRUE - use MEF-scaled axes
			%		params  (struct)		Biexp axes parameters
			%								(see Transforms.lin2logicle)
            %
            %   If inputs are not given, default to normal biexponential for 
            %   both X and Y axes with axis labels on.
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
            			
            % Check inputs
			zCheckInputs_biexpAxes();
			
            % Set axes limits
			if doMEF
				
				AXES_MAX = 2^18 * Transforms.MEF_CONVERSION_FACTOR;
				AXES_MIN = -1.5e2 * Transforms.MEF_CONVERSION_FACTOR;
				
				% Tick values
				minorTickVals = Transforms.lin2logicle( sort( ...
					[-10^5, -(1:9).*10^4, 0, ...
					 (1:9).*10^4, (1:9).*10^5, ...
					 (1:9).*10^6, (1:9).*10^7, (1:2).*10^8]), ...
					 doMEF, params);

				% Tick labels (took out -10^6 in first position)
				text = { '-10^5', '', '', '', '', '', '', '', '', '', '', ...
						 '', '', '', '', '', '', '', '', '', '10^5', ...
						 '', '', '', '', '', '', '', '', '10^6', ...
						 '', '', '', '', '', '', '', '', '10^7', ...
						 '', '', '', '', '', '', '', '', '10^8', ''};
					 
			else
				
				AXES_MAX = 2^18;
				AXES_MIN = -1.5e2;
				
				% Tick values
				minorTickVals = Transforms.lin2logicle( sort( ...
					[-10^2, -(1:9).*10^1, 0, ...
					 (1:9).*10^1, (1:9).*10^2, ...
					 (1:9).*10^3, (1:9).*10^4, (1:2).*1e5]), ...
					 doMEF, params);

				% Tick labels
				text = { '-10^2', '', '', '', '', '', '', '', '', '', '', ...
						 '', '', '', '', '', '', '', '', '', '10^2', ...
						 '', '', '', '', '', '', '', '', '10^3', ...
						 '', '', '', '', '', '', '', '', '10^4', ...
						 '', '', '', '', '', '', '', '', '10^5', ''};
			
			end	
			
			axTransformed = Transforms.lin2logicle( ...
					[AXES_MIN, AXES_MAX], doMEF, params);
			
			% Write on X-axis
			if (biexpX)
				set(ax, ...
				   'XLim', axTransformed, ...
				   'XTick', minorTickVals, ...
				   'XTickLabel', text) 
			end
            
			% Write on Y-axis
			if (biexpY)
				set(ax, ...
				   'YLim', axTransformed, ...
				   'YTick', minorTickVals, ...
				   'YTickLabel', text)
			end
            
			% Write on Z-axis
			if (biexpZ)
				set(ax, ...
				   'ZLim', axTransformed, ...
				   'ZTick', minorTickVals, ...
				   'ZTickLabel', text)
			end
			
            % Set remaining axes properties
			set(ax, ...
				'TickLength', [0.02, 0.025], ...
				'TickDir', 'out', ...
				'LineWidth', 1, ...
				'box', 'off')
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_biexpAxes()

				% Ensure boolean inputs
				biexpX = (~exist('biexpX', 'var') || all(logical(biexpX)));
				biexpY = (~exist('biexpY', 'var') || all(logical(biexpY)));
				biexpZ = (exist('biexpZ', 'var') && all(logical(biexpZ)));
				
				doMEF = (exist('doMEF', 'var') && all(logical(doMEF)));
				
				if exist('params', 'var')
					validateattributes(params, {'struct'}, {}, mfilename, 'params', 5);
				else
					params = struct();
				end
				
				% Set parameters
				params = Transforms.checkLogicleParams(doMEF, params);
			end
		end
		
		
		function ax = biexpAxes2(ax, scales, axProperties, params)
            % Tranforms the given axes to a scale used for biexponential views.
            %   
			%	ax = biexpAxes(ax, scales, params);
			%
            %   Inputs (all optional - defaults to X/Y both logicle using scale factor = 1)
            %       
			%		scales		<struct> (Optional) Property-value pairs of scaling 
			%					factors to use for each logicle-transformed dimension. 
			%					 - Valid fields: 'x', 'y', 'z' (case-insensitive)
			%					 - Any missing field will cause that dimension
			%					   to be plotted linearly
			%					 - The values here are used as the MEF value for
			%					   logicle conversions and are used as references 
			%					   for determining axes labels
			%
			%		axProperties	<struct> (optional) A struct of axes properties
			%						that can be given to overwrite default options
			%
			%		params		<struct> (Optional) Biexp transform parameters
			%					 - Note that the MEF parameter is overwritten
			%					   based on the scale provided for each dimension
			%					 - See Transforms.lin2logicle for more info
			%		
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
            			
            % Check inputs
			zCheckInputs_biexpAxes2();
			
			dims = reshape(fieldnames(scales), 1, []);
			for dim = dims
				
				% Don't set dim to @lower because its reference in 'scales'
				% needs to be preseved
				if ~any(strcmpi(dim{:}, {'x', 'X', 'y', 'Y', 'z', 'Z'})), continue, end
				
				% Set major ticks based on closest round log value
				scaleFactor = 10^round(log10(scales.(dim{:})));
				exps = (1:5) + log10(scaleFactor);
				expsStr = strrep(num2str(exps), ' ', ''); % Remove white space
				
				% Set MEF parameter for logicle calculation
				dimParams = params;
				dimParams.MEF = scales.(dim{:});
				
				% Create tick values
				tickVals = Transforms.lin2logicle( sort( ...
					[-10^exps(2), -(1:9).*10^exps(1), 0, ...
					 (1:9).*10^exps(1), (1:9).*10^exps(2), ...
					 (1:9).*10^exps(3), (1:9).*10^exps(4), (1:2).*10^exps(5)]), ...
					 true, dimParams);
				
				% Create tick labels
				% --> To add back in the 0, replace the last element of row 1
				% with the 0 that is commented out, then add another '' to the
				% second row and delete the final '' in the last row
				text = {['-10^', expsStr(2)], '', '', '', '', '', '', '', '', '', '', ...
						 ...'', '', '', '', '', '', '', '', '', '', '', ... '0^{ }', ...
						 '', '', '', '', '', '', '', '', '', [' 10^', expsStr(2)], ...
						 '', '', '', '', '', '', '', '', [' 10^', expsStr(3)], ...
						 '', '', '', '', '', '', '', '', [' 10^', expsStr(4)], ...
						 '', '', '', '', '', '', '', '', [' 10^', expsStr(5)], ''};
				
				% Determine axes limits
				AXES_MAX = 2^18 * scaleFactor;
				AXES_MIN = -1.5e2 * scaleFactor;
				axTransformed = Transforms.lin2logicle( ...
						[AXES_MIN, AXES_MAX], true, dimParams);
				
				% Set axes parameters
				set(ax, [dim{:}, 'Lim'], axTransformed, ...
						[dim{:}, 'Tick'], tickVals, ...
						[dim{:}, 'TickLabel'], text)
			end
			
            % Set remaining axes properties
			set(ax, 'TickLength', [0.02, 0.025], ...
					'TickDir', 'out', ...
					'LineWidth', 1, ...
					'box', 'off')
			
			Plotting.setAxProps(ax, axProperties);
			
			
				
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_biexpAxes2()

				% Check scales input
				if exist('scales', 'var')
					validateattributes(scales, {'struct'}, {}, mfilename, 'scales', 1);
				else
					scales = struct('x', 1, 'y', 1);
				end
				
				if ~exist('axProperties', 'var') || ~isstruct(axProperties)
					axProperties = struct();
				end
				
				% Check params input
				if exist('params', 'var')
					validateattributes(params, {'struct'}, {}, mfilename, 'params', 3);
				else
					params = struct();
				end
				
				% Set parameters (default to AFU scaling - so doMEF input is false)
				params = Transforms.checkLogicleParams(false, params);
			end
		end
		
		
		function ax = logAxes(ax, limits, axProperties, options)
			% Converts the axes to log scales and overwrites the ticks
			% so that they don't disappear and only show log-multiples of 5
			%
			%	ax = logAxes(ax, limits)
			%
			%	Inputs
			%
			%		ax				<Axes> The axes object to modify
			%
			%		limits			<struct> A struct where fields are 'x', 'y',
			%						and/or 'z' and the values are log10-transformed 
			%						axes limits, eg [3, 8] for [10^3, 10^8].
			%
			%		axProperties	<struct> (optional) A struct of axes properties
			%						that can be given to overwrite default options
			%
			%		options			<char, cell> (Optional) Logical flags:
			%							'zero': Sets the bottom tick value and
			%							label to '0' - not a transformation!
			%							'loglin': Plots log ticks and such for
			%							log-transformed (but linear) data
			%
			%	Outputs
			%	
			%		ax				<Axes> The modified axes
			
			zCheckInputs_logAxes();
			
			dims = reshape(fieldnames(limits), 1, []);
			for di = 1:numel(dims)
				
				dim = dims{di};
				
				% Don't set dim to @lower because its reference in 'scales'
				% needs to be preseved
				if ~any(strcmpi(dim, {'x', 'X', 'y', 'Y', 'z', 'Z'})), continue, end
				
				% If more than 7 log decades, skip every ceil(N/7)
				drange = limits.(dim)(1) : ...
						 ceil((limits.(dim)(2) - limits.(dim)(1)) / 7) : ...
						 limits.(dim)(2);
				
				ticklabels = cellfun(@(x) ['10^{' num2str(x), '}'], ...
							num2cell(drange), 'uniformoutput', false);
				if ismember('zero', options)
					ticklabels(1) = {'0'};
				end
				
				if ismember('loglin', options)
					set(ax, [dim, 'scale'], 'linear', ...
							[dim, 'lim'], limits.(dim), ...
							[dim, 'tick'], drange, ...
							[dim, 'ticklabel'], ticklabels);
				else
					set(ax, [dim, 'scale'], 'log', ...
							[dim, 'lim'], 10.^limits.(dim), ...
							[dim, 'tick'], 10.^(drange), ...
							[dim, 'ticklabel'], ticklabels)
				end
			end
			
			% Set remaining axes properties
			set(ax, 'TickLength', [0.02, 0.025], ...
					'TickDir', 'out', ...
					'LineWidth', 1, ...
					'box', 'off')
			
			Plotting.setAxProps(ax, axProperties);
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_logAxes()
				validateattributes(ax, {'matlab.graphics.axis.Axes'}, {}, mfilename, 'ax', 1);
				validateattributes(limits, {'struct'}, {}, mfilename, 'limits', 2);
				
				if ~exist('axProperties', 'var') || ~isstruct(axProperties)
					axProperties = struct();
				end
				
				if exist('options', 'var')
					validateattributes(options, {'char', 'cell'}, {}, mfilename, 'options', 4);
					if ischar(options), options = {options}; end % for convenience
					options = cellfun(@lower, options, 'uniformoutput', false);
				else
					options = {};
				end
			end
						
		end
		
		
		function cbar = biexpColorbar(ax, doMEF, params)
			% Generates a colorbar w/ biexponential labels for a given axes
			%
			%	cbar = biexpColorbar(ax)
			%
			%	Optional inputs:
			%
			%		doMEF	(logical)		TRUE - use MEF-scaled axes
			%		params  (struct)		Biexp axes parameters
			%								(see Transforms.lin2logicle)
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
			
			zCheckInputs_biexpColorbar()
			
			cbar = colorbar('peer', ax, 'EastOutside');
			
			% Set axes limits
% 			cbar.LimitsMode = 'manual';
% 			AXES_MAX = 2^18;
% 			AXES_MIN = -1.5e2;
% 			cbar.Limits = Transforms.lin2logicle([AXES_MIN, AXES_MAX]);
		    
			if doMEF
				% Tick values
				tickVals = Transforms.lin2logicle( sort( ...
					[-10^6, -(1:9).*10^5, -(1:9).*10^4, 0, ...
					 (1:9).*10^4, (1:9).*10^5, (1:9).*10^6, ...
					 (1:9).*10^7, (1:9).*10^8, 1e9]), doMEF, params);
				
				% Tick labels
				text = {'-10^6', '', '', '', '', '', '', '', '', '', ...
						 '', '', '', '', '', '', '', '', '', '  0^{ }', '', ...
						 '', '', '', '', '', '', '', '', '', ...
						 '', '', '', '', '', '', '', '', '10^6', ...
						 '', '', '', '', '', '', '', '', '10^7', ...
						 '', '', '', '', '', '', '', '', '10^8', ...
						 '', '', '', '', '', '', '', '', '10^9'};
			else
				% Tick values
				tickVals = Transforms.lin2logicle( sort( ...
					[-10^2, -(1:9).*10^1, 0, ...
					 (1:9).*10^1, (1:9).*10^2, ...
					 (1:9).*10^3, (1:9).*10^4, 1e5]), doMEF, params);
            
				% Tick labels
				text = {'-10^2', '', '', '', '', '', '', '', '', '', '  0^{ }', '', ...
						 '', '', '', '', '', '', '', '', '10^2', ...
						 '', '', '', '', '', '', '', '', '10^3', ...
						 '', '', '', '', '', '', '', '', '10^4', ...
						 '', '', '', '', '', '', '', '', '10^5'};
			end
			
			cbar.Ticks = tickVals;
			cbar.TickLabels = text;
			cbar.TickDirection = 'out';
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_biexpColorbar()
				
				doMEF = (exist('doMEF', 'var') && all(logical(doMEF)));
				
				if exist('params', 'var')
					validateattributes(params, {'struct'}, {}, mfilename, 'params', 5);
				else
					params = struct();
				end
				
				% Set parameters
				params = Transforms.checkLogicleParams(doMEF, params);
				
			end
		end
		
		
		function cbar = biexpColorbar2(ax, scale, params)
            % Generates a colorbar w/ biexponential labels for a given axes.
            %   
			%	cbar = biexpColorbar2(ax, scale, params);
			%
            %   Inputs 
			%
			%		ax			<axes> Axes handle to plot onto
            %       
			%		scale		<numeric> (optional) scaling factor to use for the 
			%					colorbar's logicle-transformation
			%					 - The values here are used as the MEF value for
			%					   logicle conversions and are used as references 
			%					   for determining axes labels
			%					 - Defaults to 1
			%
			%		params		<struct> (Optional) Biexp transform parameters
			%					 - Note that the MEF parameter is overwritten
			%					   based on the scale provided for each dimension
			%					 - See Transforms.lin2logicle for more info
			%
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
			
			% Check inputs
			zCheckInputs_biexpColorbar2();
			
			cbar = colorbar('peer', ax, 'EastOutside');
			
			% Adjust to closest round log value
			scaleFactor = 10^round(log10(scale));
			exps = (1:5) + log10(scaleFactor);
			expsStr = strrep(num2str(exps), ' ', ''); % Remove white space

			% Set MEF parameter for logicle calculation
			cbarParams = params;
			cbarParams.MEF = scaleFactor;

			% Create tick values
			tickVals = Transforms.lin2logicle( sort( ...
				[-10^exps(2), -(1:9).*10^exps(1), 0, ...
				 (1:9).*10^exps(1), (1:9).*10^exps(2), ...
				 (1:9).*10^exps(3), (1:9).*10^exps(4), (1:2).*10^exps(5)]), ...
				 true, cbarParams);

			% Create tick labels
			text = { ['-10^', expsStr(2)], '', '', '', '', '', '', '', '', '', '', ...
					 '', '', '', '', '', '', '', '', '', ['10^', expsStr(2)], ...
					 '', '', '', '', '', '', '', '', ['10^', expsStr(3)], ...
					 '', '', '', '', '', '', '', '', ['10^', expsStr(4)], ...
					 '', '', '', '', '', '', '', '', ['10^', expsStr(5)], ''};
			
			cbar.Ticks = tickVals;
			cbar.TickLabels = text;
			cbar.TickDirection = 'out';
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_biexpColorbar2()

				% Check scales input
				if exist('scale', 'var')
					validateattributes(scale, {'numeric'}, {'positive'}, mfilename, 'scale', 1);
				else
					scale = 1;
				end
				
				% Check params input
				if exist('params', 'var')
					validateattributes(params, {'struct'}, {}, mfilename, 'params', 2);
				else
					params = struct();
				end
				
				% Set parameters (default to AFU scaling - so doMEF input is false)
				params = Transforms.checkLogicleParams(false, params);
			end
		end

        
		function [numElements, binCenters] = biexhist(Y, numBins, showPlot)
			%BIEXHIST Logiclly-scaled (Parks,et al.) histogram. 
			%   BIEXHIST(Y) bins the elements of Y into 25 equally spaced containers
			%   and produces a histogram plot of the results.  If Y is a
			%   matrix, hist works down the columns.
			% 
			%   N = hist(Y, numBins, showPlot), where numBins is a scalar, uses numBins bins.
			%	
			%	showPlot is an optional logical flag to plot the histogram.
			%
			%   Example:
			%       [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs('sample.fcs');
			%       GreenChannel = getChannel(fcshdr,'FIT');
			%       GreenData = fcsdat(:,GreenChannel);
			%       biexhist(greenData)
			%
			% Written By
			% Breanna DiAndreth
			% bstillo@mit.edu
			% Weiss Lab, MIT
			%
			% Update Log:
			%   2015-10-31:		Changed ylimits to reflect FlowJo.  
			%					Changed default number of bins
			
			if (~exist('numBins', 'var')), numBins = 25; end
			
			ylog = Transforms.lin2logicle(Y);
			sh = range(ylog) / 20;
			binEdges = linspace(min(ylog) - sh, max(ylog) + sh, numBins + 1);
			binCenters = zeros(1, numel(binEdges) - 1);
			for b = 1:(numel(binCenters))
				binCenters(b) = mean([binEdges(b), binEdges(b + 1)]);
			end
			
			numElements = histcounts(ylog, binEdges);
			
%             % Normalize to max
%             maxN = round(max(nelements) * 1e-3) / 1e-3;
% 			  if maxN < max(nelements)
% 				  maxN = maxN + 1000;
%             end            
%             nelements=nelements./max(nelements).*100;

			if (exist('showPlot', 'var') && showPlot)
				area(binCenters, numElements, 'FaceColor', [0.5 0.5 0.5], 'LineWidth', 1.5);
				Plotting.biexpAxes(gca(), true, false);
			end
        end
        
        
        function handle = histFit(data, channel, dataType, faceColor, axPosition)
            % Plots the data as a histogram alongside the probability density functions
            %
            % - optional input dataType sets the type of data (raw, comp), used for plotting
            %
            % - optional input faceColor sets the color of the histogram faces
            %
            % - optional input axPosition sets spacing between histograms. Useful to increase size of
            % plots for best visuals.
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
            
            % Check existence of optional input dataType, assign default
            %   faceColor and axPosition are also optional, but they are looked for later
            if (~exist('dataType', 'var'))
                dataType = 'raw';
            end
            
            % Check inputs
            zCheckInputs_histFit(data, channel, dataType);
            
            % Find data sizes
            [H, W, D] = size(data);
            xrange = linspace(-2, 5, 1e5);
            if isfield(data(1).(channel), 'mus')
                P = numel(data(1).(channel).mus);
                shapes = 'ox+*sdv^<>ph';
                assert(P < length(shapes), 'Number of means must be fewer than %d\n', length(shapes));
            end
            
            % Do plotting
            handle = figure();
            spIdx = 0;
            for i = 1:H
                for j = 1:W
                    
                    % Setup subplots
                    spIdx = spIdx + 1;
                    ax = subplot(H, W, spIdx);
                    hold(ax, 'on')
                    
                    for k = 1:D
                        d = data(i, j, k).(channel).(dataType);
                        d = d(d > 1e-2);
                        h = histogram(ax, log10(d), 'Normalization', 'pdf');
                        
                        % Modify axes
                        ax.XScale = 'linear';
                        ax.YTickLabel = {};
                        if (exist('faceColor', 'var'))
                            h.FaceColor = faceColor;
                        end
%                         if (exist('axPosition', 'var'))
%                             ax.Position([3, 4]) = axPosition;
%                         end
                        ax.XLim = [-2, 5];
                        
                        % Add pdf and mean markers if fit previously
                        if isfield(data(i, j, k).(channel), 'pdf')
                            plot(xrange, data(i, j, k).(channel).pdf, 'linewidth', 2)
                        end
                        if isfield(data(i, j, k).(channel), 'mus')
                            for p = 1:P
                                plot(data(i, j, k).(channel).mus(p), 0, strcat('k', shapes(p)))
                            end
                        end
                    end
                end
            end
            
            
            % --- Helper Function --- %
            
            
            function zCheckInputs_histFit(data, channel, dataType)
                % Checks the inputs to make sure they are valid
                validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                validatestring(channel, fieldnames(data(1)), mfilename, 'channel', 2);
                validatestring(dataType, fieldnames(data(1).(channel)), mfilename, 'dataType', 3);
            end
        end
		
        
        function scatterBarPlot(ax, data, channel, dataType, colorMaps, transformation)
            % Creates a bar chart with data plotted in scatter format - gives exact representation
            % of the data essentially. 
            %
            %   Inputs:
            %       ax              The axes to plot on
            %       data            A standard struct with the data to be plotted
            %                        - Must have 'channel' and dataType given by user
            %                        - elements are plotted in a row in order
            %       channel         The channel to plot data from
            %                        - Note: The x-axis distribution is randomly normally selected
            %       dataType        The data type to use, eg 'raw', or 'scComp'
            %       colorMaps       (optional) A cell array of ColorMap objects - determines 
            %                       the colormap used for density representation
            %                        - Defaults to parula if no input given
            %       transformation  (optional) The type of transformation to make of the data
            %                        - options: 'log10', 'logicle' (default)
            %
            %
            % Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
            
            % Check inputs
            zCheckInputs_scatterBarPlot(data, channel, dataType);
            
            % Ensure axis is held on
            hold(ax, 'on');
            
            % Set colorbar for density plotting
			if (~exist('colorMaps', 'var'))
				colorMaps = cell(size(data));
				colorMaps(:) = {ColorMap('parula')};
			end
            
            % Plot data
            for i = 1:numel(data)
                                
                % Create a random normal distribution to plot points on the X-axis
                
                % Assign x- and y-axis datasets
                if exist('transformation', 'var')
                    switch transformation
                        case 'log10'
                            ydata = log10(data(i).(channel).(dataType));
%                             ydata = ydata(ydata > 4.25);
                        case 'logicle'
                            ydata = Transforms.lin2logicle(data(i).(channel).(dataType));
                    end
                else
                    ydata = Transforms.lin2logicle(data(i).(channel).(dataType));
                end
                
                valid = ~(isnan(ydata) | isinf(ydata));
                ydata = ydata(valid);
                xrnd = randn(numel(ydata), 1);
                xdata = i + (xrnd - mean(xrnd)) / range(xrnd);
                
                % Plot points
                Plotting.densityplot(ax, xdata, ydata, 10000, 'hist', colorMaps{i});
            end
            
            % Set Y-axis to biexp
            if ~exist('transformation', 'var') || strcmpi(transformation, 'logicle')
                Plotting.biexpAxes(ax, false, true);
            end
            ax.XTick = 1:numel(data);
            ax.XTickLabelRotation = -45;
            ax.FontSize = 14;
            
            
            % --- Helper Function --- %
            
            
            function zCheckInputs_scatterBarPlot(data, channel, dataType)
                % Checks the inputs to make sure they are valid
                validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                validatestring(channel, fieldnames(data(1)), mfilename, 'channel', 2);
                validatestring(dataType, fieldnames(data(1).(channel)), mfilename, 'dataType', 3);
            end
        end
        
        
        function batchBoxPlot(ax, data, channel, dataType, numPoints, transformation)
            % Creates a box plot for the data. 
            %   Automatically converts rows in data to separate groupings
            %
            %   Inputs:
            %       ax              The axes to plot on
            %       data            A standard struct with the data to be plotted
            %                        - Must have 'channel' and dataType given by user
            %                        - elements are plotted in a row in order
            %       channel         The channel to plot data from
            %                        - Note: The x-axis distribution is randomly normally selected
            %       dataType        The data type to use, eg 'raw', or 'scComp'
            %       numPoints       (optional) The # of points to use for making the box 
            %                        - Default: 6,000
            %       transformation  (optional) The type of transformation to make of the data
            %                        - options: 'log10' (default), 'logicle'
            %
            %
            % Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
            
            % Check inputs
            zCheckInputs_batchBoxPlot(data, channel, dataType);
            
            % Check number of points to use for plotting, default to 6,000 if no number given
            if exist('numPoints', 'var')
                validateattributes(numPoints, {'numeric'}, {'integer'}, mfilename, 'numPoints', 5);
            else
                numPoints = 6000;
            end
            
            % Ensure axis is held on
            hold(ax, 'on');
            
            % Plot data
            dataMatrix = zeros(numPoints, numel(data));
            for i = 1:numel(data)
                
                % Assign y-axis dataset
                ydata = data(i).(channel).(dataType);
                if exist('transformation', 'var')
                    switch transformation
                        case 'log10'
                            ydata = log10(ydata(ydata > 0));
                        case 'logicle'
                            ydata = Transforms.lin2logicle(ydata);
                    end
                else
                    ydata = Transforms.lin2logicle(data(i).(channel).(dataType));
                end
                
                valid = ~(isnan(ydata) | isinf(ydata));
                ydata = ydata(valid);
                
                % Ensure there are enough points.
                assert(numel(ydata) >= numPoints, ...
                    'Must have at least %d points!\nOnly have %d in data set %d', ...
                    numPoints, numel(ydata), i);
                
                % Extract numPoints random points from ydata into the data matrix
                dataMatrix(:, i) = ydata(randperm(numel(ydata), numPoints));
            end
            
            % Make a compact graph if more than 10 bars
            if numel(data) <= 10
                plotStyle = 'traditional';
            else
                plotStyle = 'compact';
            end
            
%             % Determine groupings by rows in data
%             numGroups = size(data, 1);
%             groupings = repmat(1:numGroups, [1, size(data, 2)]);
            
            % Make boxplot
            whos dataMatrix
            boxplot(dataMatrix, 'plotstyle', plotStyle)
            
            % Set Y-axis to biexp if applicable
            if ~exist('transformation', 'var') || strcmpi(transformation, 'logicle')
                Plotting.biexpAxes(ax, false, true);
            end
            ax.XTick = 1:numel(data);
            ax.XTickLabelRotation = -45;
            ax.FontSize = 14;
            
            
            % --- Helper Function --- %
            
            
            function zCheckInputs_batchBoxPlot(data, channel, dataType)
                % Checks the inputs to make sure they are valid
                validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                validatestring(channel, fieldnames(data(1)), mfilename, 'channel', 2);
                validatestring(dataType, fieldnames(data(1).(channel)), mfilename, 'dataType', 3);
            end
        end

		
        function batchViolinPlot(ax, data, channel, dataType, colors)
            % Creates a standard violin plot for each element in the given data struct, plotted
            % in sequence.
            %
            %   Inputs:
            %       ax              The axes to plot on
            %       data            A standard struct with the data to be plotted
            %                        - Must have 'channel' and dataType given by user
            %                        - elements are plotted in a row in order
            %       channel         The channel to plot data from
            %       dataType        The data type to use, eg 'raw', or 'scComp'
            %       colors          (optional) An Nx3 array of rbg color triplets - determines 
            %                       the color order used for the face of each violin
            %                        - N must be numel(data)
            %                        - Defaults to standard MATLAB sequence if no input given
            %
            % Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
            
            % Check inputs
            zCheckInputs_batchViolinPlot(data, channel, dataType);
            
            % Ensure axis is held on
            hold(ax, 'on');
            
            % Set colorbar for density plotting
            if (exist('colors', 'var'))
                assert(size(colors, 1) == numel(data));
            else
                colors = repmat([ ...
                    0       0.4470  0.7410
                    0.8500  0.3250  0.0980
                    0.9290  0.6940  0.1250
                    0.4940  0.1840  0.5560
                    0.4660  0.6740  0.1880
                    0.3010  0.7450  0.9330
                    0.6350  0.0780  0.1840], ...
                    ceil(numel(data) / 7), 1);
            end
            
            % Plot data
			binEdges = 0:0.1:4.5;
            for di = 1:numel(data)
                
                % Assign x- and y-axis datasets
                ydata = Transforms.lin2logicle(data(di).(channel).(dataType));
                
                % Plot points
                Plotting.violinplot(ax, ydata, di, binEdges, colors(di, :));
            end
            
            % Set Y-axis to biexp
            Plotting.biexpAxes(ax, false, true);
            ax.XTick = 1:numel(data);
            ax.XTickLabelRotation = -45;
            ax.FontSize = 14;
            
            
            % --- Helper Function --- %
            
            
            function zCheckInputs_batchViolinPlot(data, channel, dataType)
                % Checks the inputs to make sure they are valid
                validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                validatestring(channel, fieldnames(data(1)), mfilename, 'channel', 2);
                validatestring(dataType, fieldnames(data(1).(channel)), mfilename, 'dataType', 3);
            end
        end
        
        
        function ax = lineDensityPlot(inputs)
            % Creates a line density plot for each element in the given data struct. Lines are
            % essentially the connected quantities of bins from a histogram of a channel's data. 
            %
            %   Inputs (struct array w/ fields below):
            %       ax          (optional) The axes to plot on - default is gca()
            %       data        A standard struct with the data to be plotted
            %                    - Must have 'channel' and 'dataType' given by user
            %                    - elements are plotted in a row in order
            %       channel     The channel to plot data from
            %       dataType    The data type to use, eg 'raw', or 'scComp'
            %       colors      (optional) An Nx3 array of rbg color triplets - determines 
            %                   the color order used for the face of each violin
            %                    - N must be numel(data)
            %                    - Defaults to standard MATLAB sequence if no input given
            %       xscale      (optional) String: 'linear', 'log', 'logicle' 'logicleMEF'
			%					Determines the scale of the x-axis. Default = 'logicle'
			%		yscale		(optional) String: 'log', 'linear'
			%					Determines the scale of the y-axis. Default = 'log'
			%		options		<cell, char> (optional) A cell array of strings (or a single 
			%					string) specifying optional plotting behavior:
			%						'shade'		Flag to shade in the area under the line. 
			%						'counts'	Flag to plot bin counts rather than PDF
			%
			%	Outputs
			%		ax				The axes handle for the generated plot
            %
            % Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
            
            % Check inputs
            [data, channel, dataType] = zCheckInputs_lineDensityPlot(inputs);
            
            % Ensure axis is held on
            if isfield(inputs, 'ax')
                ax = inputs.ax;
            else
                ax = gca();
            end
            hold(ax, 'on');
            
            % Set line color order
			if (isfield(inputs, 'colors'))
				colors = inputs.colors;
				assert(size(colors, 1) == numel(data), ...
					   'Size of colors incorrect');
			else
				colors = repmat(lines(7), ceil(numel(data) / 7), 1);
			end
            
			if isfield(inputs, 'xscale')
				validatestring(inputs.xscale, ...
					{'linear', 'lin', 'log', 'logicle', 'biexp', 'logicleMEF', 'MEF'}, ...
					mfilename, 'options.xscale');
			else
				inputs.xscale = 'logicle'; % Default behavior
			end
			
            % Plot data
            for i = 1:numel(data)
                
				% Extract data
				xdata = data(i).(channel).(dataType);
				switch inputs.xscale
					case {'linear', 'lin'}
						break % no change necessary
					case {'log'}
						xdata = log10(xdata);
					case {'logicle', 'biexp'}
						xdata = Transforms.lin2logicle(xdata);
					case {'logicleMEF', 'MEF'}
						xdata = Transforms.lin2logicle(xdata, true);
					otherwise
						error('Unrecognized xscale given: %s', inputs.xscale);
				end
				
				% Plot data
				if (isfield(inputs, 'shade') && inputs.shade)
					options = {'shade'};
				else
					options = {};
				end
				Plotting.singleLineDensity(ax, xdata, colors(i, :), options)
				
            end
            
            % Scale axes
			switch inputs.xscale
				case {'linear', 'lin'}
					ax.XScale = 'linear';
				case {'log'}
					ax.XScale = 'log';
				case {'logicle', 'biexp'}
					Plotting.biexpAxes(ax, true, false);
				case {'logicleMEF', 'MEF'}
					Plotting.biexpAxes(ax, true, false, false, true);
			end
			
			% Set Y-axis
			if isfield(options, 'yscale')
				validatestring(options.yscale, {'linear', 'log'}, mfilename, 'options.yscale')
			else
				options.yscale = 'log';
			end
            ax.YScale = options.yscale;
            
			ax.FontSize = 16;
             
            
            % --- Helper Function --- %
            
            
            function [data, channel, dataType] = zCheckInputs_lineDensityPlot(inputs)
                
                % Ensure necessary inputs are present
                reqFields = {'data', 'channel', 'dataType'};
				missingFields = setdiff(reqFields, fieldnames(inputs));
                assert(isempty(missingFields), ...
						'Missing field: %s', missingFields{:})
                
                % Extract data
                data = inputs.data;
                channel = inputs.channel;
                dataType = inputs.dataType;
                
                % Checks the inputs to make sure they are valid
                validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                validatestring(channel, fieldnames(data(1)), mfilename, 'channel', 2);
                validatestring(dataType, fieldnames(data(1).(channel)), mfilename, 'dataType', 3);
            end
		end
        
		
		function [h, binCenters, binCounts] = singleLineDensity(ax, xdata, color, options, edges, float)
			% Plots a sinlge line representing a histogram of values in a vector
			% of given values. The Y-axis represents the probability density function.
			%
			%	Plotting.singleLineDensity(ax, xdata, color, options)
			%
			%	Inputs
			%		ax			<axes> axes handle to plot on
			%		xdata		<numeric> A vector of data values to plot
			%		color		<numeric> (optional) A single 1x3 RGB color 
			%					triplet to color the line 
			%		options		<cell, char> (optional) A cell array of strings (or a single 
			%					string) specifying optional plotting behavior:
			%						'shade'		Flag to shade in the area under the line. 
			%						'counts'	Flag to plot bin counts rather than probability
			%						'maxnorm'	Flag to normalize by the max count/probability 
			%									value (has precedence).
			%						'kde'		Flag to use kernel density estimation, not histogram
			%						'smooth'	Flag to smooth the histogram data
			%						'boxplot'	Flag to add lines denoting percentiles: 
			%										[5, 25, 50 (median), 75, 95]
			%		edges		<numeric> (optional) Edges to use for histogram binning
			%					Default: 30 evenly spaced bins chosen by histcounts()
			%		float		<numeric> (optional) A value to add to the baseline 
			%					of the histogram to 'float' it above the x-axis
			%
			%	Outputs
			%		h			<handle> The handle to the plotted element (eg for legends)
			%		binCenters  <numeric> An array of computed bin center points
			%		binCounts   <numeric> An array of computed bin counts
			%						(normalized and smoothed if requested in inputs)
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
			
			% TODO
			%	1) Rename to histogram
			%	2) Check if bins are log or linearly spaced to compute centers
			%	and maybe to also pre-transform the data
			
			
			zCheckInputs_singleLineDensity();
			
			% Fix negative infinite values by setting the resulting values to the 
			% otherwise minimum value.
			% Fix NaN and inf values by removing them.
			if any(isnan(xdata) | isinf(xdata))
				warning('Inf/NaN values detected - removing')
				valid = ~(isnan(xdata) | isinf(xdata));
				xdata = xdata(valid);
			end
			
			% Automatically find number of bins based on data
			NUM_BINS = 32; % Power of 2 for KDE
			if ismember('kde', options)
				% Use kernel to estimate density
				[~, binCounts, binCenters] = kde(xdata, NUM_BINS);
			else
				% Use histogram to estimate density
				if ismember('counts', options)
					histNormMode = 'count';
				else
					histNormMode = 'pdf';
				end
				if exist('edges', 'var')
					binCounts = histcounts(xdata, edges, ...
							'normalization', histNormMode);
				else
					[binCounts, edges] = histcounts(xdata, NUM_BINS, ...
							'normalization', histNormMode);
				end
				if ismember('maxnorm', options)
					binCounts = binCounts ./ max(binCounts) * 0.8;
				end
				
				% Estimate bin centers by averaging the edges
				binCenters = zeros(numel(edges) - 1, 1);
				for j = 2:length(edges)
					binCenters(j - 1) = mean(edges([j - 1, j]));
				end
			end
			
			% Ignore bins w/o data, otherwise gaps in the line plot will show up
			binCounts = reshape(binCounts, [], 1);
			hasPoints = (binCounts > 0);
			binCounts = binCounts(hasPoints);
			binCenters = binCenters(hasPoints);
			if ismember('smooth', options)
				binCounts = smooth(binCounts);
			end
			
			% Adjust lines to a common baseline
			minX = min(binCenters);
			maxX = max(binCenters);
			minY = min(binCounts) / 10 + float; % In case log axes, don't use 0
			binCenters = [minX; binCenters; maxX];
			binCounts = [minY; binCounts + float; minY];
			
			% Plot w/ or w/o shading
			if ismember('shade', options)
				patch(ax, binCenters, binCounts, color, 'edgeAlpha', 0,  'facealpha', 0.4);
			end
			if ismember('boxplot', options)
				pcts = prctile(xdata, [5, 25, 50, 75, 95]);
				linewidths = [0.5, 1, 2, 1, 0.5];
				for pi = 1:numel(pcts)
					[~, ci] = min(abs(binCenters - pcts(pi)));
					plot(ax, binCenters([ci, ci]), [float, binCounts(ci)], ...
							'color', color, 'linewidth', linewidths(pi))
				end
			end
			h = plot(ax, binCenters, binCounts, 'color', color, 'linewidth', 2);
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_singleLineDensity()
				
				validateattributes(ax, {'matlab.graphics.axis.Axes'}, {}, mfilename, 'ax', 1);
				validateattributes(xdata, {'numeric'}, {'vector'}, mfilename, 'xdata', 2);
				
				% If no color is given, use the next default color
				if ~exist('color', 'var')
					color = ax.ColorOrder(ax.ColorOrderIndex, :);
				end
				validateattributes(color, {'numeric'}, {'vector'}, mfilename, 'color', 3);
				
				% Color should be 1x3 or 3x1
				assert(length(color) == 3, 'Input ''color'' must be length 3')
				color = reshape(color, 1, []);
				
				% If no optionsa are given, make an empty input set
				if ~exist('options', 'var') || isempty(options)
					options = {};
				end
				validateattributes(options, {'cell', 'char'}, {}, mfilename, 'options', 4);
				if ischar(options), options = {options}; end % For simplicity
				options = cellfun(@lower, options, 'uniformoutput', false);
				
				if exist('edges', 'var')
					validateattributes(edges, {'numeric'}, {}, mfilename, 'edges', 5);
				end
				
				if exist('float', 'var')
					validateattributes(float, {'numeric'}, {'scalar'}, mfilename, 'float', 6);
				else
					float = 0;
				end
			end
		end
		
		
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
			%							errorbars. To change the calculations,
			%							use the 'options' input. 
            %
 			%		options			<cell, char> (Optional) A cell array of strings (or 
			%						a single string) specifying optional plotting behavior:
			%							'geomean'	Compute bar height w/ geomean 
			%										rather than mean.
			%							'geostd'	Compute errorbar height w/
			%										geostd rather than std.
			%							'relerr'	Compute errorbars using
			%										relative error (for log plots)
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
				% --> Mode of sorted xpos to handle cases where xpos is not for
				%	  all samples that will eventually be plotted and where the
				%	  xpos are out of order
				if (numel(xpos) > 1)
					dx = mode(diff(sort(xpos))); 
				else
					dx = 1;
				end
				
				xvals = linspace(0.7, 1.3, size(yvals, 1))';
				rxi = zeros(size(yvals));
				for xi = 1:numel(xpos)
					rxi(:, xi) = randperm(numel(xvals))';
				end
				
				xvals = repmat(xpos, size(yvals, 1), 1) + dx .* (xvals(rxi) - 1);
				
				scatter(ax, xvals(:), yvals(:), 10, 'k', 'filled')
			end
				
			Plotting.setAxProps(ax, axProperties);
			
			
			% --- Helper Functions --- %
			
			
			function [ymeans, ystdsNeg, ystdsPos, errProperties] = zCheckInputs_barPlot()
				
				% Force into row vector
				xpos = reshape(xpos, 1, []);
				assert(size(xpos, 2) == size(yvals, 2), ...
						'Number of columns in ''xpos'' and ''yvals'' inputs must match');
				
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
				
				if ismember('relerr', options)
					[ystdsNeg, ystdsPos] = relError(ymeans, ystds);
				else
					ystdsNeg = ystds;
					ystdsPos = ystds;
				end
				
				if ismember ('unierr', options)
					% Make errors only show up on the appropriate 'side' of the bar
					ystdsNeg = ystdsNeg .* (ymeans < 0);
					ystdsPos = ystdsPos .* (ymeans > 0);
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
        
        function [ax, h, cbar] = standardHeatmap(dataMatrix, rowLabels, colLabels, cmap, options)
            % Generates a standard heatmap with fontsize 14 and data arranged nicely
            %
            %   Inputs:
            %
            %       dataMatrix      The data to be plotted (matrix). Normally, heatmaps flip
            %                       the rows for some reason (like an image) - this is 
            %                       corrected such that the heatmap "starts" in the top left
            %
            %       rowLabels       Cell list of strings with row labels
			%						(the function accounts for flipping as mentioned above). 
            %
            %       colLabels       Cell list of strings with col labels 		
			%
			%		cmap			<numeric, char, ColorMap> (Optional) 
			%						The colormap to represent data values. 
			%						- Can input an Nx3 matrix of RGB values, a
			%						   ColorMap object pre-initialized with a
			%						   color, or a string indicating which
			%						   ColorMap to initialize. 
			%						 - ColorMap and char inputs will yield 100
			%						   unique color values on the given scale
			%						 - Defualt = parula(100);
			%
			%		options			(optional) Struct of additional optional inputs
			%			norm            The dimension to normalize on
            %								1 or 'column' = within columns
            %								2 or 'row' = within rows
            %								3 or 'none' = no normalization
			%			dataType		The type of data being plotted
			%								'mef'/'mefl' will plot onto standard biexpMEF colorbar
			%								'raw' will plot onto standard bixep colorbar
			%								'null' (defualt) will not adjust the colobar
			%			position		A 1x4 position description for the figure to be plotted
			%			doCluster       Boolean value which flags to create a clustergram
            %							rather than a heatmap. Clustergram is created with standard
            %							euclidean distance calculation for both rows and cols.
            %			symmetric       Boolean value which flags the heatmap to be 
            %							symmetric around zero. (Default = False)
			%			colpdist		The distance metric to use for columns
			%								[If using dendrograms]
			%			rowpdist		The distance metric to use for rows
			%								[If using dendrograms]
			%
            %   Outputs: 
            %
            %       ax              A handle to the figure axes created
            %       h               A handle to the heatmap/clustergram object
            % 
            % Written By 
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
            %
            % Update Log:
            %   2016-05-02:		Added symmetric argument
			%	2018-01-15:		Changed extra inputs to options, added colorbar editing
			%
			
			% Check optional inputs exist
			if ~exist('options', 'var')
				options = struct();
			end
			
            % Invert data rows because stupid heatmap
            flippedData = flipud(dataMatrix);
            
            % Create heatmap/clustergram
			if (isfield(options, 'doCluster') && options.doCluster)
				doCluster = true;
				h = clustergram(flippedData);
				if isfield(options, 'rowpdist')
					h.RowPDist = options.rowpdist;
				else
					h.RowPDist = 'correlation';
				end
				if isfield(options, 'colpdist')
					h.ColumnPDist = options.colpdist;
				else
					h.ColumnPDist = 'correlation';
				end
			else
				doCluster = false;
				h = HeatMap(flippedData);
			end
			
			if isfield(options, 'norm')
				h.Standardize = upper(options.norm);
				if isfield(options, 'range')
					h.DisplayRange = options.range; 
				end
			end
			
			if isfield(options, 'symmetric')
				h.Symmetric = options.symmetric;
			else
				h.Symmetric = false;
			end
			
			if exist('cmap', 'var')
				h.Colormap = Plotting.checkCmap(cmap);
			else
				h.Colormap = parula(100);
			end
			
            h.RowLabels = fliplr(rowLabels);
            h.ColumnLabels = colLabels;
            
%             h.view();
            
            % Create figure/axes object and update font size / add colorbar
            ax = h.plot();
			if (~doCluster)
				% The colorbar throws off the dendrogram, so it needs to be generated in other ways
				if ismember('dataType', fieldnames(options))
					switch options.dataType
						case {'raw'}
							cbar = Plotting.biexpColorbar(ax);
						case {'mef', 'mefl'}
							cbar = Plotting.biexpColorbar(ax, true);
						otherwise
							cbar = colorbar('peer', ax, 'EastOutside');
							cbar.TickDirection = 'out';
							% No major changes to colorbar, just move ticks out
					end
				end
			end
			
			% Change figures size/position
			if ismember('position', fieldnames(options))
				set(gcf(), 'position', options.position);
			end
			
            ax.FontSize = 14;
		end
		
		
		function ax = heatmap(ax, dataMatrix, edges, labels, cmap, axProperties, options)
			% Creates a heatmap representing the given data.
			%
			%	ax = heatmap(ax, dataMatrix, edges, labels, cmap, axProperties, options)
			%
			%	 - 2D data is plotted as a simple heatmap on a single axis
			%	 - 3D data is plotted as a series of 2D heatmaps stacked in
			%	   the 3rd dimension on a single axis
			%
			%	Inputs
			%
			%		ax				<Axes> Handle to the axes to plot onto
			%
			%		data			<numeric> A 2-3 dimensional matrix containing
			%						data to be plotted.
			%
			%		edges			<cell> A cell list of numeric arrays containing 
			%						the edge values for each heatmap element.
			%						Each element of the cell list corresponds
			%						with one dimension. 
			%
			%		labels			<cell> A cell list of labels for each dimension.
			%						NOTE: The last entry is the label for the
			%						colorbar. Thus, # labels = ndims(data) + 1
			%
			%		cmap			<numeric, char, ColorMap> (Optional) 
			%						The colormap to represent data values. 
			%						 - Can input an Nx3 matrix of RGB values, a
			%						   ColorMap object pre-initialized with a
			%						   color, or a string indicating which
			%						   ColorMap to initialize. 
			%						 - ColorMap and char inputs will yield 100
			%						   unique color values on the given scale
			%						 - Defualt = viridis(100);
			%
			%		axProperties	<struct> (Optional) Property-value pairs for
			%						the axes properties
			%
			%		options			<struct> Optional property-value pairs:
			%							'biexp', <dimensions> enables biexponential 
			%							  axes in the given dimensions 
			%								({'C', 'X', 'Y', 'Z'} accepted, where 
			%								'C' corresponds w/ the colorbar)
			%							  - Can also pass a struct with the
			%								dimension names as fields to set the
			%								MEF scaling factor for axes labels
			%							'doMEF' If TRUE, does logicle conversion
			%							  with MEF-unit scaling (1e3)
			%								(default = FALSE)
			%							'params': <params> enables setting the
			%							  logicle function parameters
			%							  (see Transforms.lin2logicle())
			%								(defaults to 'standard' values
			%							'min': <min val> enables setting the
			%							  lower bound for color-data conversion
			%								(default = 0)
			%							'max': <max val> enables setting the
			%							  upper bound for color-data conversion
			%								(default = 4.5)
			%							'fig': <fig handle> enables setting the
			%							  figure to plot onto
			%								(default = new fig)
			%							'categorical': <logical> A flag to plot
			%							  the axes 'categorically', ie with no
			%							  tick lines, 90deg rotated x-axis labels,
			%							  and larger font size
			%								(default = FALSE)
			%							'plotColorbar': <logical> A flag to plot
			%							  the colorbar.
			%								(default = TRUE)
			%
			%	Outputs
			%
			%		figBinHmap		<handle> A handle to the generated figure 
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%			
			
			doTitle = xCheckInputs_heatmap();
			
			for d3 = 1:size(dataMatrix, 3)
				
				% Setup patch coordinates
				patchesX = zeros(4, size(dataMatrix, 1) * size(dataMatrix, 2));
				patchesY = zeros(size(patchesX));
				patchesZ = ones(size(patchesX)) * (d3 - 1);
				patIdx = 0;

				% Adjust for differential edges in higher dims. Linear indexing
				% ensures that this works no matter how the edges were formatted
				if iscell(edges{2})
					currSubsIdx = {d3, options.d4, options.d5};
					currLinearIdx = sub2ind(options.dataSize(3:end), currSubsIdx{:});
					e2 = edges{2}{currLinearIdx};
				else
					e2 = edges{2};
				end

				% Extract patch coordinates
				for d2 = 1:size(dataMatrix, 2)

					% Adjust first dimension for dims 2+
					if iscell(edges{1})
						currSubsIdx = {d2, d3, options.d4, options.d5};
						currLinearIdx = sub2ind(options.dataSize(2:end), currSubsIdx{:});
						e1 = edges{1}{currLinearIdx};
					else
						e1 = edges{1};
					end

					for d1 = 1:size(dataMatrix, 1)

						patIdx = patIdx + 1;

						patchesX(:, patIdx) = [
							e1(d1); e1(d1 + 1); e1(d1 + 1); e1(d1)];

						patchesY(:, patIdx) = [
							e2(d2); e2(d2); e2(d2 + 1); e2(d2 + 1)];

					end
				end

				% Get colors for each patch
				colors = Plotting.getColors(dataMatrix(:, :, d3), cmap, options);

				% Plot all patches
				patch(ax, patchesX, patchesY, patchesZ, ...
					reshape(colors, [size(patchesX, 2), 1, 3]), ...
					'EdgeColor', 'none', 'FaceColor', 'flat');

				% Set axes properties/labels
				% - We do this before setting values on the Z-axis since
				%	we want any 3rd dim bins to be evenly spaced
% 						Plotting.biexpAxes(ax, ismember('X', options.biexp), ...
% 							   ismember('Y', options.biexp), ...
% 							   ismember('Z', options.biexp), ...
% 							   options.doMEF, options.params);
				Plotting.logAxes(ax, options.log, [], options.logParams);
				Plotting.biexpAxes2(ax, options.biexp, [], options.params);
				xlabel(ax, labels{1}, 'fontsize', 10);
				ylabel(ax, labels{2}, 'fontsize', 10);
			end
			
			% Plot Z (3D) labels if applicable
			if (size(dataMatrix, 3) > 1 && doTitle(3))
				centers3 = cell(1, size(dataMatrix, 3));
				for c3i = 1:numel(centers3)
					if iscell(edges{3}(c3i))
						centers3{c3i} = edges{3}{c3i};
					else % Numeric edges
						if (size(dataMatrix, 3) < numel(edges{3}))
							edgeVal = mean([edges{3}(c3i), edges{3}(c3i + 1)]);
						else % Centers given as input
							edgeVal = edges{3}(c3i);
						end
						if edges{3}(c3i) > 1e3
							% Show large numbers in scientific notation
							centers3{c3i} = sprintf('%.2g', edgeVal);
						else
							centers3{c3i} = num2str(edgeVal);
						end
					end
				end
				
				% Need to fix ticks, otherwise they will be 
				% added/removed as plot size changes
				ax.ZTick = 0:(size(dataMatrix, 3) - 1); 
				ax.ZTickLabel = centers3;
				ax.ZLim = [0, size(dataMatrix, 3) - 1];
				zlabel(ax, labels{3});
			end
			
			% Standard axes properties
			if options.categorical
				set(ax, 'XTickLabelRotation', 90, ...	% Vertical X labels
						'TickLength', [0, 0], ...		% No ticks
						'fontsize', 10, ...				% Readable font size
						'XTick', 0.5:size(dataMatrix, 1), ... % Center the
						'YTick', 0.5:size(dataMatrix, 2), ... % category labels
						'XLim', [0, size(dataMatrix, 1)], ... % Sometimes it auto-adjusts wrong
						'YLim', [0, size(dataMatrix, 2)]);
			end
			
			% Settable properties
			Plotting.setAxProps(ax, axProperties);
			
			% Colormap
			if options.plotColorbar
				fO = fieldnames(options.biexp);
				[~, idxO] = ismember({'c', 'C'}, fO);
				if any(idxO)
	% 				cbar = Plotting.biexpColorbar(ax, options.doMEF, options.params);
					fC = fO{idxO(find(idxO > 0, 1))};
					cbar = Plotting.biexpColorbar2(ax, options.biexp.(fC), options.params);
				else
					cbar = colorbar(ax);
	% 				cbar.Limits = [min(dataMatrix(:)), max(dataMatrix(:))];
				end
				colormap(ax, cmap);
				cbar.Limits = [options.min, options.max];
				ax.CLim = cbar.Limits; % Need this or the color won't fill the entire bar
				cbar.Label.String = labels{end};
				if isfield(axProperties, 'FontSize')
					cbar.Label.FontSize = axProperties.FontSize;
				else
					cbar.Label.FontSize = 10;
				end
			end
			
			% --- Helper Functions --- %
			
			
			function doTitle = xCheckInputs_heatmap()
				
				% Check axes
				validateattributes(ax, {'matlab.graphics.axis.Axes'}, {}, mfilename, 'ax', 1);
				
				% Check data + dimensions
				validateattributes(dataMatrix, {'numeric'}, {}, mfilename, 'data', 2);
				assert(ndims(dataMatrix) > 1, 'Data must have at least two dimensions!')
				assert(ndims(dataMatrix) < 4, 'Data must have no more than 3 dimensions!')
				
				% Determine size of data
				dataSize = ones(1, 3); % 5 Max allowed number of dimensions
				dataSizeTemp = size(dataMatrix);
				dataSize(1:numel(dataSizeTemp)) = dataSizeTemp;
				
				% Check edges
				if isempty(edges)
					edges = cell(1, ndims(dataMatrix));
					for di = 1:ndims(dataMatrix)
						if di <= 2
							edges{di} = 0:size(dataMatrix, di);
						else
							edges{di} = [];
						end
					end
				else
					validateattributes(edges, {'cell'}, {}, mfilename, 'edges', 3);
				end
				assert(numel(edges) == ndims(dataMatrix), ...
						'Number of set of edges (%d) does not match dimensions of data! (%d)', ...
						numel(edges), ndims(dataMatrix));
				doTitle = true(size(edges));
				for ei = 1:numel(edges)
					if isempty(edges{ei})
						doTitle(ei) = false;
						continue
					end
					if ei <= 2
						if iscell(edges{ei}) % Cell list of edges for each higher dim
							assert(numel(edges{ei}) == prod(dataSize(ei + 1:end)), ...
									'Number of dim %d edges variants does not match data', ei);
							for dei = 1:numel(edges{ei})
								assert(numel(edges{ei}{dei}) == (size(dataMatrix, ei) + 1), ...
									'Number of edges supplied for dim %d is incorrect!', ei)
							end
						else
							assert(numel(edges{ei}) == (size(dataMatrix, ei) + 1), ...
								   'Number of edges supplied for dim %d is incorrect!', ei)
						end
					else % Allow specifying exact positions for dims 3+
						assert(numel(edges{ei}) == size(dataMatrix, ei) || ...
							   numel(edges{ei}) == (size(dataMatrix, ei) + 1), ...
							   'Number of edges supplied for dim %d is incorrect!', ei)
					end
				end
				
				% Check labels
				validateattributes(labels, {'cell'}, {}, mfilename, 'labels', 4);
				assert(numel(labels) == (ndims(dataMatrix) + 1), ...
					'Number of labels (%d) does not match data dimensionality + 1! (%d)', ...
					numel(labels), ndims(dataMatrix));
				
				% Check colormap
				if exist('cmap', 'var')
					validateattributes(cmap, {'numeric', 'char', 'ColorMap'}, {}, mfilename, 'cmap', 5);
					if ischar(cmap)
						cmap = ColorMap(cmap).getColormap(100);
					elseif isa(cmap, 'ColorMap')
						cmap = cmap.getColorMap(100);
					end
				else
					cmap = viridis(100);
				end
				
				% Check axes properties
				if (exist('axProperties', 'var') && ~isempty(axProperties))
					validateattributes(axProperties, {'struct'}, {}, mfilename, 'axProperties', 6);
				else
					axProperties = struct();
				end
				
				% Check options
				if exist('options', 'var')
					validateattributes(options, {'struct'}, {}, mfilename, 'options', 7);
				else
					options = struct();
				end
				if ~isfield(options, 'log'), options.log = {}; end
				if ~isfield(options, 'logParams'), options.logParams = {}; end
				if ~isfield(options, 'biexp'), options.biexp = {}; end
				if ~isfield(options, 'params'), options.params = struct(); end
				if ~isfield(options, 'doMEF'), options.doMEF = false; end
				if options.doMEF, autoscale = 1e3; else, autoscale = 1; end
				if ~isfield(options, 'min'), options.min = 0; end
				if ~isfield(options, 'max'), options.max = 4.5; end
				if ~isfield(options, 'categorical'), options.categorical = false; end
				if ~isfield(options, 'plotColorbar'), options.plotColorbar = true; end
				
				% Options added by Plotting.binHeatmap()
				if ~isfield(options, 'd4'), options.d4 = 1; end
				if ~isfield(options, 'd5'), options.d5 = 1; end
				if ~isfield(options, 'dataSize'), options.dataSize = dataSize; end
				
				% Convert biexp/log options to struct if necessary
				if ischar(options.biexp), options.biexp = {options.biexp}; end
				if iscell(options.biexp)
					biexp = struct();
					for i = 1:numel(options.biexp)
						try
							biexp.(options.biexp{i}) = autoscale;
						catch % in case the fieldname is invalid
						end
					end
					options.biexp = biexp;
				end
				if iscell(options.log)
					logopt = struct();
					for i = 2:numel(options.log)
						try % Assume a cell list of axes-limits were passed
							logopt.(options.log{i - 1}) = options.log{i};
						catch % in case the fieldname is invalid or otherwise fails
						end
					end
					options.log = logopt;
				end
			end
		end
		
		
		function figBinHmap = binHeatmap(dataMatrix, edges, labels, cmap, axProperties, options)
			% Creates a heatmap representing the given data. The presentation 
			% depends on the dimensionality of the data.
			%
			%	figBinHmap = binHeatmap(dataMatrix, edges, labels, cmap, axProperties, options)
			%
			%	 - 2D data is plotted as a simple heatmap on a single axis
			%	 - 3D data is plotted as a series of 2D heatmaps stacked in
			%	   the 3rd dimension on a single axis
			%		* To plot individual Z-stacks on different axes, pass the
			%		  3rd dimension as a singleton and push the data to 4/5. 
			%	 - 4D data os plotted as a series of 3D data in multiple
			%	   subplots occupying one row of the figure. 
			%	 - 5D data is plotted like 4D data, with the 5th dimension
			%	   occupying individual subplot rows of the figure. 
			%
			%	Inputs
			%
			%		data			<numeric> A 2-5 dimensional matrix containing
			%						data to be plotted.
			%
			%		edges			<cell> A cell list of numeric arrays containing
			%						the edge values for each heatmap element.
			%						Each element of the cell list corresponds
			%						with one dimension. 
			%						 - For dimensions 3+, tick values can be
			%						   given instead of edges
			%						 - Edges for dimensions 3+ can also be
			%						   categorical, given as a cell list of
			%						   strings with length equal to size(data, N)
			%						 - If the contents of one array element is
			%						   itself a cell list (which must be the
			%						   same size as dimensions 3+), then each
			%						   element will be used independently for
			%						   each higher-dimensional heatmap. 
			%
			%		labels			<cell> A cell list of labels for each dimension.
			%						NOTE: The last entry is the label for the
			%						colorbar. Thus, # labels = ndims(data) + 1
			%
			%		cmap			<numeric, char, ColorMap> (Optional) 
			%						The colormap to represent data values. 
			%						 - Can input an Nx3 matrix of RGB values, a
			%						   ColorMap object pre-initialized with a
			%						   color, or a string indicating which
			%						   ColorMap to initialize. 
			%						 - ColorMap and char inputs will yield 100
			%						   unique color values on the given scale
			%						 - Defualt = viridis(100);
			%
			%		axProperties	<struct> (Optional) Property-value pairs for
			%						the axes properties for each axes object
			%						plotted by the function. 
			%
			%		options			<struct> Optional property-value pairs:
			%							'biexp', <dimensions> enables biexponential 
			%							  axes in the given dimensions 
			%								({'C', 'X', 'Y', 'Z'} accepted, where 
			%								'C' corresponds w/ the colorbar)
			%							  - Can also pass a struct with the
			%								dimension names as fields to set the
			%								MEF scaling factor for axes labels
			%							'doMEF' If TRUE, does logicle conversion
			%							  with MEF-unit scaling
			%								(default = FALSE)
			%							'params': <params> enables setting the
			%							  logicle function parameters
			%							  (see Transforms.lin2logicle())
			%								(defaults to 'standard' values)
			%							'min': <min val> enables setting the
			%							  lower bound for color-data conversion
			%								(default = 0)
			%							'max': <max val> enables setting the
			%							  upper bound for color-data conversion
			%								(default = 4.5)
			%							'fig', <fig handle> enables setting the
			%							  figure to plot onto
			%								(defualt = new fig)
			%							'categorical': <logical> A flag to plot
			%							  the axes 'categorically', ie with no
			%							  tick lines, 90deg rotated x-axis labels,
			%							  and larger font size
			%								(default = FALSE)
			%							'plotColorbar': <logical> A flag to plot
			%							  the colorbar.
			%								(default = TRUE)
			%
			%	Outputs
			%
			%		figBinHmap		<Figure> A handle to the generated figure 
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
			%	2018-03-06		Added option for giving figure handle to plot on
			
			% TODO: 
			%	1) Make it so N or N+1 edges can be used for D1-2 as well 
			%		(automatically treat as categorical if so)
			%	2) Rename to batchHeatmap or something like that
			
			[figBinHmap, doTitle] = zCheckInputs_binHeatmap();
			
			spIdx = 0;
			for d5 = 1:size(dataMatrix, 5)
				for d4 = 1:size(dataMatrix, 4)
					
					spIdx = spIdx + 1;
					ax = subplot(size(dataMatrix, 5), size(dataMatrix, 4), spIdx);
					
					if (size(dataMatrix, 3) == 1)
						elSubs = [1, 2];
					else
						elSubs = [1, 2, 3];
					end
					
					optionsHM = options;
					optionsHM.d4 = d4; optionsHM.d5 = d5;
					optionsHM.plotColorbar = false;
					ax = Plotting.heatmap(ax, dataMatrix(:, :, :, d4, d5), ...
							edges(elSubs), labels([elSubs, numel(labels)]), ...
							cmap, axProperties, optionsHM);
					
					% Plot titles w/ 4/5D labels if applicable
					titleTxt = {};
					if (size(dataMatrix, 4) > 1 && doTitle(4))
						if iscell(edges{4}(d4))
							center4 = edges{4}{d4};
						else
							if (size(dataMatrix, 4) < numel(edges{4}))
								edgeVal = mean([edges{4}(d4), edges{4}(d4 + 1)]);
							else
								edgeVal = edges{4}(d4);
							end
							if edges{4}(d4) > 1e3
								% Show large numbers in scientific notation
								center4 = sprintf('%.2g', edgeVal);
							else
								center4 = num2str(edgeVal);
							end
						end
						titleTxt = [titleTxt; {sprintf('%s = %s', labels{4}, center4)}]; %#ok<AGROW>
					end
					if (size(dataMatrix, 5) > 1 && doTitle(5))
						if iscell(edges{5}(d5))
							center5 = edges{5}{d5};
						else
							if (size(dataMatrix, 5) < numel(edges{5}))
								edgeVal = mean([edges{5}(d5), edges{5}(d5 + 1)]);
							else
								edgeVal = edges{5}(d5);
							end
							if edges{5}(d5) > 1e3
								% Show large numbers in scientific notation
								center5 = sprintf('%.2g', edgeVal);
							else
								center5 = num2str(edgeVal);
							end
						end
						titleTxt = [titleTxt; {sprintf('%s = %s', labels{5}, center5)}]; %#ok<AGROW>
					end
					title(ax, titleTxt);
					
% 					% Standard axes properties
% 					if options.categorical
% 						set(ax, 'XTickLabelRotation', 90, ...	% Vertical X labels
% 								'TickLength', [0, 0], ...		% No ticks	
% 								'fontsize', 10, ...				% Readable font size
% 								'XTick', 0.5:size(dataMatrix, 1), ... % Center the
% 								'YTick', 0.5:size(dataMatrix, 2), ... % category labels
% 								'XLim', [0, size(dataMatrix, 1)], ... % Sometimes it auto-adjusts wrong
% 								'YLim', [0, size(dataMatrix, 2)]);
% 					end
% 					
% 					% Settable properties
% 					Plotting.setAxProps(ax, axProperties);
				end
			end
			
			if options.plotColorbar
				fO = fieldnames(options.biexp);
				[~, idxO] = ismember({'c', 'C'}, fO);
				if any(idxO)
	% 				cbar = Plotting.biexpColorbar(ax, options.doMEF, options.params);
					fC = fO{idxO(find(idxO > 0, 1))};
					cbar = Plotting.biexpColorbar2(ax, options.biexp.(fC), options.params);
				else
					cbar = colorbar(ax);
	% 				cbar.Limits = [min(dataMatrix(:)), max(dataMatrix(:))];
				end
				colormap(ax, cmap);
				cbar.Limits = [options.min, options.max];
				cbar.FontSize = ax.FontSize;
				ax.CLim = cbar.Limits; % Need this or the color won't fill the entire bar
				cbar.Label.String = labels{end};
				cbar.Label.FontSize = ax.XLabel.FontSize;
				if isfield(axProperties, 'FontSize')
					cbar.Label.FontSize = axProperties.FontSize;
				else
					cbar.Label.FontSize = 10;
				end
			end
			
			
			% --- Helper Functions --- %
			
			
			function [figBinHmap, doTitle] = zCheckInputs_binHeatmap()
				
				% Check data + dimensions
				validateattributes(dataMatrix, {'numeric'}, {}, mfilename, 'data', 1);
				assert(ndims(dataMatrix) > 1, 'Data must have at least two dimensions!')
				assert(ndims(dataMatrix) < 6, 'Data must have no more than 5 dimensions!')
				
				% Determine size of data
				dataSize = ones(1, 5); % 5 Max allowed number of dimensions
				dataSizeTemp = size(dataMatrix);
				dataSize(1:numel(dataSizeTemp)) = dataSizeTemp;
				
				% Check edges
				if isempty(edges)
					edges = cell(1, ndims(dataMatrix));
					for di = 1:ndims(dataMatrix)
						if di <= 2
							edges{di} = 0:size(dataMatrix, di);
						else
							edges{di} = [];
						end
					end
				else
					validateattributes(edges, {'cell'}, {}, mfilename, 'edges', 2);
				end
				assert(numel(edges) == ndims(dataMatrix), ...
						'Number of set of edges (%d) does not match dimensions of data! (%d)', ...
						numel(edges), ndims(dataMatrix));
				doTitle = true(size(edges));
				for ei = 1:numel(edges)
					if isempty(edges{ei})
						doTitle(ei) = false;
						continue
					end
					if ei <= 2
						if iscell(edges{ei}) % Cell list of edges for each higher dim
							assert(numel(edges{ei}) == prod(dataSize(ei + 1:end)), ...
									'Number of dim %d edges variants does not match data', ei);
							for dei = 1:numel(edges{ei})
								assert(numel(edges{ei}{dei}) == (size(dataMatrix, ei) + 1), ...
									'Number of edges supplied for dim %d is incorrect!', ei)
							end
						else
							assert(numel(edges{ei}) == (size(dataMatrix, ei) + 1), ...
								   'Number of edges supplied for dim %d is incorrect!', ei)
						end
					else % Allow specifying exact positions for dims 3+
						assert(numel(edges{ei}) == size(dataMatrix, ei) || ...
							   numel(edges{ei}) == (size(dataMatrix, ei) + 1), ...
							   'Number of edges supplied for dim %d is incorrect!', ei)
					end
				end
				
				% Check labels
				validateattributes(labels, {'cell'}, {}, mfilename, 'labels', 3);
				assert(numel(labels) == (ndims(dataMatrix) + 1), ...
					'Number of labels (%d) does not match data dimensionality + 1! (%d)', ...
					numel(labels), ndims(dataMatrix));
				
				% Check colormap
				if exist('cmap', 'var')
					validateattributes(cmap, {'numeric', 'char', 'ColorMap'}, ...
							{}, mfilename, 'cmap', 4);
					if ischar(cmap)
						cmap = ColorMap(cmap).getColormap(100);
					elseif isa(cmap, 'ColorMap')
						cmap = cmap.getColorMap(100);
					end
				else
					cmap = viridis(100);
				end
				
				% Check axes properties
				if (exist('axProperties', 'var') && ~isempty(axProperties))
					validateattributes(axProperties, {'struct'}, {}, mfilename, 'axProperties', 5);
				else
					axProperties = struct();
				end
				
				% Check options
				if exist('options', 'var')
					validateattributes(options, {'struct'}, {}, mfilename, 'options', 7);
				else
					options = struct();
				end
				if ~isfield(options, 'log'), options.log = {}; end
				if ~isfield(options, 'logParams'), options.logParams = {}; end
				if ~isfield(options, 'biexp'), options.biexp = {}; end
				if ~isfield(options, 'params'), options.params = struct(); end
				if ~isfield(options, 'doMEF'), options.doMEF = false; end
				if options.doMEF, autoscale = 1e3; else, autoscale = 1; end
				if ~isfield(options, 'min'), options.min = 0; end
				if ~isfield(options, 'max'), options.max = 4.5; end
				if ~isfield(options, 'categorical'), options.categorical = false; end
				if ~isfield(options, 'plotColorbar'), options.plotColorbar = true; end
				if ~isfield(options, 'fig')
					figBinHmap = figure();
				else
					figBinHmap = options.fig;
				end
				options.dataSize = dataSize;
				
				% Convert biexp/log options to struct if necessary
				if ischar(options.biexp), options.biexp = {options.biexp}; end
				if iscell(options.biexp)
					biexp = struct();
					for i = 1:numel(options.biexp)
						try
							biexp.(options.biexp{i}) = autoscale;
						catch % in case the fieldname is invalid
						end
					end
					options.biexp = biexp;
				end
				if iscell(options.log)
					logopt = struct();
					for i = 2:numel(options.log)
						try % Assume a cell list of axes-limits were passed
							logopt.(options.log{i - 1}) = options.log{i};
						catch % in case the fieldname is invalid or otherwise fails
						end
					end
					options.log = logopt;
				end
			end
		end
		
		
		function cmap = checkCmap(cmap)
			% Simple function to check if a given colormap is an Nx3 matrix of
			% RGB colors, a colormap name, or a ColorMap object, returning an
			% Nx3 matrix of RGB colors if the input is valid. For colormap name
			% or ColorMap inputs, the generated colormap contains 100 colors.
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%
			
			validateattributes(cmap, {'numeric', 'char', 'ColorMap'}, ...
					{}, mfilename, 'cmap', 2);
			
			if ischar(cmap)
				cmap = ColorMap(cmap).getColormap(100);
			elseif strcmpi(class(cmap), 'ColorMap')
				cmap = cmap.getColormap(100);
			end
			
			if (min(cmap(:)) < 0)
				warning('Colormap min value less than 0! Shifting up...')
				cmap = cmap + min(cmap(:));
			end
			if (max(cmap(:)) > 1)
				warning('Colormap max value greater than 1! Scaling down...')
				cmap = cmap ./ max(cmap(:));
			end
		end
		
	end
	
end