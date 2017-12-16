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
    % Written/Compiled by Ross Jones
    % Weiss Lab, MIT
    % Last updated 2016-05-27
    
	methods (Static)
		
		
		function density = computeDensity(dataMatrix, mode, numPoints, nonZero)
			% Computes the density of a given set of points
			%
			%	density = computeDensity(dataMatrix, mode, numPoints, nonZero)
			%
			%	Inputs
			%		dataMatrix		<numeric> An NxM matrix of N elements in M 
			%						dimensions for which to compute density
            %
            %       mode			<char> Determines how density is calculated
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
			
			% Check inputs
			dataValid = checkInputs_computeDensity();
            numDims = size(dataMatrix, 2);
			
			% Check mode of point coloration
			subsample = randperm(size(dataValid, 1), numPoints);
			dataSub = dataValid(subsample, :);
			switch mode 
				case {'neighbors', 'fastn'}
					% Find how many neighbours there are less than dX and dY away.
					if (strcmpi(mode, 'fastn'))
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
						sh = (max(dataValid(:, d)) - min(dataValid(:, d))) / 20;
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
			
			
			function dataValid = checkInputs_computeDensity()
				
				validateattributes(dataMatrix, {'numeric'}, {}, mfilename, 'data', 1);
				validatestring(mode, {'neighbors', 'fastn', 'kernel', 'hist'}, mfilename, 'mode', 2);
				
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
		
		
		function [sortIdx, colors] = getColors(inputData, colorMap, options)
			% Generates a color representing the value of each element of a
			% given input vector using a given ColorMap object
			%
			%	[sortedInput, colors] = getColors(input, colorMap, options)
			%
			%	Inputs
			%		inputData	<numeric> A vector for which to generate colors.
			%					NaN values are set to the lowest possible value.
			%		
			%		colorMap	<ColorMap> A ColorMap object initialized with
			%					the desired colormap to grab colors from
			%
			%		options		<struct> A number of optional inputs:
			%						'numColors': An integer (default = 100)
			%						'min': The minimum color threshold 
			%							(default = min(inputData))
			%						'max': The maximum color threshold 
			%							(default = max(inputData))
			%
			%	Outputs
			%		sortIdx		The sorting order for inputData
			%
			%		colors		An Nx3 matrix of color values corresponding
			%					with each N elements in sort(inputData)
			
			% Check inputs
			[numColors, MIN, MAX] = checkInputs_getColors();
			
			% Sort so brightest cells are plotted on top
			[sortedInput, sortIdx] = sort(inputData);
			
			% Convert density to color
			cm = colorMap.getColormap(numColors);
			% The -1/+1 ensures that no value is exactly 0 or > max index
			colorIdxs = max(min(ceil((sortedInput - MIN) ./ (MAX - MIN) .* numColors), numColors), 1);
			try
				colors = cm(colorIdxs, :);
			catch ME
				MIN
				MAX
				numColors
				min(colorIdxs)
				max(colorIdxs)
				error('Encountered error getting colors')
			end
			
			
			% --- Helper Functions --- %
			
			
			function [numColors, MIN, MAX] = checkInputs_getColors()
				
				validateattributes(inputData, {'numeric'}, {'vector'}, mfilename, 'inputData', 1);
				validateattributes(colorMap, {'ColorMap'}, {}, mfilename, 'colorMap', 2);
				
				% Fix for inf/nan values
				inputData(inputData == inf) = max(inputData(~isinf(inputData)));
				inputData(inputData == -inf) = min(inputData(~isinf(inputData)));
				inputData(isnan(inputData)) = min(inputData);
				
				% Default options
				numColors = 100;
				MIN = min(inputData);
				MAX = max(inputData);
				
				if exist('options', 'var')
					
					if isfield(options, 'numColors')
						numColors = round(options.numColors);
						validateattributes(numColors, {'numeric'}, {'scalar', 'positive'}, mfilename, 'options.numColors', 3);
					end
					
					if isfield(options, 'min')
						MIN = options.min;
						validateattributes(MIN, {'numeric'}, {'scalar'}, mfilename, 'options.min', 3)
					end
					
					if isfield(options, 'max')
						MAX = options.max;
						validateattributes(MAX, {'numeric'}, {'scalar'}, mfilename, 'options.max', 3)
					end
				end
			end
		end
		
		
		function densityplot(ax, xdata, ydata, nPoints, mode, colorMap, nonZero)
			%DENSITYPLOT(XDATA,YDATA) plots the vector Y vs vector X in dot-plot form with colors
            %   of the dots indicating density
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
            %       mode (string)       Determines how density is calculated
            %                            'normal'        Computes density directly from the points
            %                            'fast'          Uses 1/10 points to compute density
            %                            'kernel'        Kernel density estimation 
            %                                            (auto-selects bandwidth - see kde2d.m)
            %                            'hist'          Interpolates density from a 2D histogram
            %                            <numerical>     Colors the points based on a given set of values
            %                                            Must be the same size as xdata/ydata
            %       
            %       colorMap (ColorMap) The ColorMap opject to use for plotting
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
			%   Written by
			%   Breanna Stillo
			%   bstillo@mit.edu
			%   Last Updated: 2014-10-14
			%
			%   Edit 2015-02-06 (Ross) - Doubled speed by slightly changing how neighbours is calculated
            %   Edit 2016-03-28 (Ross) - Added kernel density estimation and merged with fastdensity
            %                            calculation function
            
            % Check inputs
            assert(all(size(xdata) == size(ydata)), 'xdata and ydata are different sizes!')
            validateattributes(xdata, {'numeric'}, {'vector'}, mfilename, 'xdata', 2);
            validateattributes(ydata, {'numeric'}, {'vector'}, mfilename, 'ydata', 3);
            validateattributes(nPoints, {'numeric'}, {'scalar'}, mfilename, 'nPoints', 4);
            validateattributes(mode, {'char', 'numeric'}, {'vector'}, mfilename, 'mode', 5);
            if ischar(mode)
                validatestring(mode, {'normal', 'fast', 'kernel', 'hist'}, mfilename, 'mode', 5);
            else
                assert(all(size(mode) == size(xdata)), 'If numeric, mode must be the same size as the data!')
            end
            validateattributes(colorMap, {'ColorMap'}, {}, mfilename, 'colorMap', 6);
			if exist('nonZero', 'var')
				nonZero = any(logical(nonZero(:)));
				validateattributes(nonZero, {'logical'}, {}, mfilename, 'nonZero', 7);
			else
				nonZero = false;
			end
			nPoints = round(nPoints);

            % Remove complex values
            x = real(xdata);
            y = real(ydata);
            
            % Fix NaN and inf values by removing them.
			if any(isnan(x) | isinf(x) | isnan(y) | isinf(y))
				warning('Inf/NaN values detected - removing')
				valid = ~(isnan(x) | isinf(x) | isnan(y) | isinf(y));
			else
				valid = true(size(x));
			end
			if (nonZero), valid = (valid & (x >= 0) & (y >= 0)); end
            
			% Check mode of point coloration
            if ischar(mode)
				
				% Only need to check x and y for valid points
				x = x(valid);
				y = y(valid);
				numPoints = min(numel(x), nPoints);
				subsample = randperm(numel(x), numPoints);
				x = x(subsample);
				y = y(subsample);
				switch mode 
					case {'fast', 'normal'}
						% Find how many neighbours there are less than dX and dY away.
						if (strcmpi(mode, 'fast'))
							skip = 10;
						else
							skip = 1;
						end

						% Find the range of each vector
						xdist = max(x) - min(x);
						ydist = max(y) - min(y);

						% Find the density of each vector by dividing the number of points by the
						% range. Dividing 100 by this number defines some distance one would expect 
						% close data points to be from one another.
						xdens = numPoints / xdist;
						ydens = numPoints / ydist;
						dX = 100 / xdens;
						dY = 100 / ydens;

						neighbors = zeros(numPoints, 1);
						for j = 1:numPoints
							xval = x(j);
							yval = y(j);

							neighbors(j) = sum( ...
								(abs(x(1:skip:end) - xval) < dX) & ...
								(abs(y(1:skip:end) - yval) < dY));
						end

	%                     % Convert neighbors to log scale to get better view of data
	%                     neighbors = round(10 * Transforms.lin2logicle(neighbors));

						[sortedNeighbors, sortIdx] = sort(neighbors);

						% Convert # neighbors to color
						nColors = max(sortedNeighbors) + 1;
						cm = colorMap.getColormap(nColors);
						colors = cm(sortedNeighbors + 1, :);

					case 'kernel'
						% Reshape input vector to column format
						x = reshape(x, [], 1);
						y = reshape(y, [], 1);

						% Estimate density with kernel
						[bw, density, meshX, meshY] = kde2d([x, y]);

	%                     % Convert density to log scale to get better view of data
	%                     density = Transforms.lin2logicle(density);

	%                     fprintf(1, 'Kernal density estimation with bandwidth: %.3f\n', bw);
						interpDensity = interp2(meshX, meshY, density, x, y);

						[sortedDensity, sortIdx] = sort(interpDensity);

						% Convert density to color
						nColors = 100;
						cm = colorMap.getColormap(nColors);
						% The -1/+1 ensures that no value is exactly 0 or > max index
						MIN = min(sortedDensity);
						if (nonZero), MIN = max(0, MIN); end
						colors = cm(ceil((sortedDensity - MIN) ...
										./ (max(sortedDensity) - MIN) ...
										.* (nColors - 1) + 1), :);

					case 'hist'
						% Reshape input vector to column format
						x = reshape(x, [], 1);
						y = reshape(y, [], 1);
						nBins = 25;

						% Use histogram to calculate true density - select # bins based on # points
						edgesX = linspace(median(x) - 5 * std(x), median(x) + 5 * std(x), nBins + 1);
						edgesY = linspace(median(y) - 5 * std(y), median(y) + 5 * std(y), nBins + 1);
						[binCounts, edgesX, edgesY] = histcounts2(x, y, edgesX, edgesY);

						% Estimate bin centers by averaging the edges
						binCoordX = zeros(numel(edgesX) - 1, 1);
						binCoordY = zeros(numel(edgesY) - 1, 1);
						for i = 2:length(edgesX)
							binCoordX(i - 1) = mean(edgesX([i - 1, i]));
						end
						for i = 2:length(edgesY)
							binCoordY(i - 1) = mean(edgesY([i - 1, i]));
						end

						% Create mesh of X and Y values for the bins
						[binMeshX, binMeshY] = meshgrid(binCoordX, binCoordY);

						% Interpolate over the 2D histogram counts to make a PDF for the data points
						interpDensity = interp2(binMeshX, binMeshY, binCounts', x, y);

						% Sort so brightest cells are plotted on top
						[sortedDensity, sortIdx] = sort(interpDensity);
						validInterp = ~(isnan(sortedDensity) | isinf(sortedDensity));
						sortedDensity = sortedDensity(validInterp);
						sortIdx = sortIdx(validInterp);

						% Convert density to color
						nColors = 100;
						cm = colorMap.getColormap(nColors);
						% The -1/+1 ensures that no value is exactly 0 or > max index
						MIN = min(sortedDensity);
						if (nonZero), MIN = max(0, MIN); end
						colors = cm(ceil((sortedDensity - MIN) ...
										./ (max(sortedDensity) - MIN) ...
										.* (nColors - 1) + 1), :);

				end
			else
				z = real(reshape(mode, [], 1));
				valid = (valid & ~(isnan(z) | isinf(z)));
				if (nonZero), valid = (valid & (z >= 0)); end
				x = x(valid);
				y = y(valid);
				z = z(valid);
				
				% Reduce number of points to speed density calculation
				numPoints = min(numel(x), nPoints);
				subsample = randperm(numel(x), numPoints);
				x = x(subsample);
				y = y(subsample);
                z = z(subsample);       
                [sortedZ, sortIdx] = sort(z);
				
                nColors = 100;
                cm = colorMap.getColormap(nColors);
				colorIdx = ceil(sortedZ ./ 4.5 .* (nColors - 1) + 1);
				colorIdx = max(min(colorIdx, nColors), 0);
                colors = cm(colorIdx, :);
% 				colors = cm(sortedZ, :);
            end
            
            % Plot data in sorted order so highest density points are plotted last.
            scatter(ax, x(sortIdx), y(sortIdx), 8, colors, 'filled')
        end
        
        
        function violinplot(ax, ydata, xcenter, mode, faceColor)
            %VIOLINPLOT(DATA, MODE, COLORMAP) plots the data as a violin plot, which represents 
            % density as a vertical, horizontally symmetric histogram
			%   
            %   Inputs
            %   
			%       ydata (vector)       values of data to plot
            %
            %       xcenter (integer)   The value to center the x-values on
            %
            %       mode (string)       Determines how density is calculated
            %                            'normal'        Computes density directly from the points
            %                            'fast'          Uses 1/10 points to compute density
            %                            'kernel'        Kernel density estimation 
            %                                            (auto-selects bandwidth - see kde.m)
            %                            'hist'          Interpolates density from a 2D histogram
            %       
            %       faceColor (rgb)     The color of the violin face. The face is plotted with an
            %                           alpha value of 0.5 to look nice.
            %
			%
			%
			%   Written by
			%   Ross Jones
            %   jonesr18@mit.edu
			%   Weiss Lab, MIT
			%   Last Updated: 2016-04-04
            
            % Check inputs
            checkInputs_violinplot(ydata, xcenter, mode, faceColor); 
            xcenter = round(xcenter(1)); % Ensure integer and only one point
            
            % Remove complex values
            y = real(ydata);
            
            % Fix negative infinite values by setting the resulting values to the 
            % otherwise minimum value.
            if any(y == -inf)
                warning('Negative values detected - setting to min value')
                y(y == -inf) = min(y(y ~= -inf));
            end
            
            % Reduce number of points to speed density calculation
            numPoints = min(numel(y), 5000);
            subsample = randperm(numel(y), numPoints);
            y = y(subsample);
            
            switch mode
                case {'fast', 'normal'}
                    % Find how many neighbours there are less than dD away.
                    if (strcmpi(mode, 'fast'))
                        skip = 10;
                    else
                        skip = 1;
                    end
                    
                    % Find the range of each vector
                    dRange = max(y) - min(y);

                    % Find the density of each vector by dividing the number of points by the
                    % range. Dividing 100 by this number defines some distance one would expect 
                    % close data points to be from one another.
                    yDens = numel(y) / dRange;
                    dY = 100 / yDens;
                    
                    neighbors = zeros(numPoints, 1);
                    for j = 1:numPoints
                        neighbors(j) = sum((abs(y(1:skip:end) - y(j)) < dY));
                    end
                    
                    % Sub-sample d and neighbors to reduce number of points again for plotting a
                    % more smooth violin
                    smallNumPoints = min(numel(y));
                    smallSubsample = randperm(numel(y), smallNumPoints);
                    [points, sortIdx] = sort(y(smallSubsample));
                    density = neighbors(sortIdx);
                    
                case 'kernel'
                    % Reshape input vector to column format
                    y = reshape(y, [], 1);
                    
                    % Estimate density with kernel
                    [bw, density, points] = kde(y, 512);
                    points = reshape(points, [], 1);
                    
%                     % Convert density to log scale to get better view of data
%                     density = Transforms.lin2logicle(density);
                    
%                     fprintf(1, 'Kernal density estimation with bandwidth: %.3f\n', bw);
                
                case 'hist'
                    % Reshape input vector to column format
                    y = reshape(y, [], 1);
                    
                    % Use histogram to calculate true density
                    [count, edges] = histcounts(y);
                    density = reshape(count, [], 1);
                    
%                     % Convert counts to log scale to get better view of data
%                     binCounts = Transforms.lin2logicle(binCounts);
                    
                    % Estimate bin centers by averaging the edges
                    points = zeros(length(edges) - 1, 1);
                    for i = 2:length(edges)
                        points(i - 1) = mean(edges([i - 1, i]));
                    end
            end
            
            % Adjust density scale
            density = Transforms.lin2logicle(density);
            density = interp1([min(density), max(density)], [0, 0.4], density);
            
            % Plot violins
            axes(ax);
            fill( [xcenter + density; xcenter - flipud(density)], ...
                  [points; flipud(points)], ...
... %                   faceColor, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
                  faceColor, 'EdgeColor', 'none');
           
            % --- Helper Function --- %
            
            
            function checkInputs_violinplot(ydata, xcenter, mode, faceColor)
              
                validateattributes(ydata, {'numeric'}, {'vector'}, mfilename, 'data', 2);
                validateattributes(xcenter, {'numeric'}, {}, mfilename, 'xcenter', 3);
                validatestring(mode, {'normal', 'fast', 'kernel', 'hist'}, mfilename, 'mode', 4);
                validateattributes(faceColor, {'numeric'}, {'vector'}, mfilename, 'faceColor', 5);
                assert(length(faceColor) == 3);
            end
        end
        
        
		function biexplot(x, y, plotArgs, options)
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
			%   Written by
			%   Breanna Stillo
			%   bstillo@mit.edu
			%   Last Updated: 2014-10-14;
			%
			%   Edited 2-22-16 by Ross Jones
			%   Cleaned up some code
            
            if (~exist('ax', 'var'))
                ax = gca();
            end
            
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
		
        
        function biexpAxes(ax, biexpX, biexpY, biexpZ)
            % Tranforms the given axes to a scale used for biexponential veiws.
            %   
            %   Inputs
            %       
            %       biexpX (logical)        TRUE - make x-axis biexponential
            %       biexpY (logical)        TRUE - make y-axis biexponential
            %       biexpZ (logical)        TRUE - make z-axis biexponential
            %
            %   If inputs are not given, default to biexponential for both X 
            %   and Y axes with axis labels on.
            
% 			axes(ax)
% 			ax2 = axes;
			
            % Ensure boolean inputs
            if (exist('biexpX', 'var'))
                biexpX = logical(biexpX);
            else
                biexpX = true;
            end
            if (exist('biexpY', 'var'))
                biexpY = logical(biexpY);
            else
                biexpY = true;
            end
			if (exist('biexpZ', 'var'))
				biexpZ = logical(biexpZ);
			else
				biexpZ = false;
			end
			
            % Set axes limits
			AXES_MAX = 2^18;
			AXES_MIN = -1.5e2;
            axTransformed = Transforms.lin2logicle([AXES_MIN, AXES_MAX]);

            % Major tick values
%             majorTickVals = Transforms.lin2logicle([-1e2 -1e1 0 1e1 1e2 1e3 1e4 1e5]);
		    
			% Minor tick values
			minorTickVals = Transforms.lin2logicle( sort( ...
                [-10^2, -(1:9).*10^1, 0, ...
                 (1:9).*10^1, (1:9).*10^2, ...
				 (1:9).*10^3, (1:9).*10^4, 1e5]));
            
                       % Major tick labels
			
% 			text = {'-10^6', '', '  0^{ }', '', '10^6', '10^7', '10^8', '10^9'};
			text = {'-10^2', '', '', '', '', '', '', '', '', '', '  0^{ }', '', ...
					 '', '', '', '', '', '', '', '', '10^2', ...
					 '', '', '', '', '', '', '', '', '10^3', ...
					 '', '', '', '', '', '', '', '', '10^4', ...
					 '', '', '', '', '', '', '', '', '10^5'};
            
			% Write on X-axis
			if (biexpX)
				set(ax, ...
				   'XLim', axTransformed, ...
				   'XTick', minorTickVals, ...
				   'XTickLabel', text) 
% 			   set(ax2, ...
% 				   'XLim', axTransformed, ...
% 				   'XTick', minorTickVals, ...
% 				   'XTickLabel', {})
			end
            
			% Write on Y-axis
			if (biexpY)
				set(ax, ...
				   'YLim', axTransformed, ...
				   'YTick', minorTickVals, ...
				   'YTickLabel', text)
% 			   set(ax2, ...
% 				   'YLim', axTransformed, ...
% 				   'YTick', minorTickVals, ...
% 				   'YTickLabel', {})
			end
            
			% Write on Z-axis
			if (biexpZ)
				set(ax, ...
				   'ZLim', axTransformed, ...
				   'ZTick', minorTickVals, ...
				   'ZTickLabel', text)
% 			   set(ax2, ...
% 				   'YLim', axTransformed, ...
% 				   'YTick', minorTickVals, ...
% 				   'YTickLabel', {})
			end
			
            % Set remaining axes properties
% 			set(ax, ...
% 				'TickLength', [0.04 0.05], ...
% 				'TickDir', 'out', ...
% 				'LineWidth', 1.5, ...
% 				'Box', 'off', ...
% 				'FontSize', 8)
			set(ax, ...
				'TickLength', [0.02, 0.025], ...
				'TickDir', 'out', ...
				'LineWidth', 1, ...
				'box', 'off')
			
			% Make sure ax is current before returning
% 			axes(ax);
			
%             hha3=axes;
%             set(hha3,'color','none', ...
%                 'XTick',[],'XTickLabel',{},...
%                 'YTick',[],'YTickLabel',{},...
%                 'xlim', axTransformed,...

		end
		
		
		function biexpAxesMEF(ax, biexpX, biexpY, biexpZ)
            % Tranforms the given MEF unit axes to a scale used for biexponential veiws.
            %   
            %   Inputs
            %       
            %       biexpX (logical)        TRUE - make x-axis biexponential
            %       biexpY (logical)        TRUE - make y-axis biexponential
            %       biexpZ (logical)        TRUE - make z-axis biexponential
            %
            %   If inputs are not given, default to biexponential for both X 
            %   and Y axes with axis labels on.
            
			% Make sure current figure is brought forward, which is necessary
			% for adding axis layers for minor tick vals. We have to create them
			% manually because the built-in minor tick vals are auto-generated!
% 			axes(ax)
% 			ax2 = axes;
			
            % Ensure boolean inputs
            if (exist('biexpX', 'var'))
                biexpX = logical(biexpX);
            else
                biexpX = true;
            end
			if (exist('biexpY', 'var'))
				biexpY = logical(biexpY);
			else
				biexpY = true;
			end
			if (exist('biexpZ', 'var'))
                biexpZ = logical(biexpZ);
            else
                biexpZ = false;
			end
            
            % Set axes limits
			AXES_MAX = 1.1e9;
			AXES_MIN = -4e5;
            axTransformed = Transforms.lin2logicleMEF([AXES_MIN, AXES_MAX]);

            % Major tick values
% 			majorTickVals = Transforms.lin2logicleMEF([-1e6 -1e5 0 1e5 1e6 1e7 1e8 1e9]);
		    
			% Minor tick values
			minorTickVals = Transforms.lin2logicleMEF( sort( ...
                [-10^6, -(1:9).*10^5, -(1:9).*10^4, 0, ...
                 (1:9).*10^4, (1:9).*10^5, (1:9).*10^6, ...
				 (1:9).*10^7, (1:9).*10^8, 1e9]));
            
            % Major tick labels
% 			text = {'-10^6', '', '  0^{ }', '', '10^6', '10^7', '10^8', '10^9'};
			text = {'-10^6', '', '', '', '', '', '', '', '', '', ...
					 '', '', '', '', '', '', '', '', '', '  0^{ }', '', ...
					 '', '', '', '', '', '', '', '', '', ...
					 '', '', '', '', '', '', '', '', '10^6', ...
					 '', '', '', '', '', '', '', '', '10^7', ...
					 '', '', '', '', '', '', '', '', '10^8', ...
					 '', '', '', '', '', '', '', '', '10^9'};
            
			% Write on X-axis
			if (biexpX)
				set(ax, ...
				   'XLim', axTransformed, ...
				   'XTick', minorTickVals, ...
				   'XTickLabel', text) 
% 			   set(ax2, ...
% 				   'XLim', axTransformed, ...
% 				   'XTick', minorTickVals, ...
% 				   'XTickLabel', {})
			end
            
			% Write on Y-axis
			if (biexpY)
				set(ax, ...
				   'YLim', axTransformed, ...
				   'YTick', minorTickVals, ...
				   'YTickLabel', text)
% 			   set(ax2, ...
% 				   'YLim', axTransformed, ...
% 				   'YTick', minorTickVals, ...
% 				   'YTickLabel', {})
			end
			
			% Write on Z-axis
			if (biexpZ)
				set(ax, ...
				   'ZLim', axTransformed, ...
				   'ZTick', minorTickVals, ...
				   'ZTickLabel', text)
% 			   set(ax2, ...
% 				   'YLim', axTransformed, ...
% 				   'YTick', minorTickVals, ...
% 				   'YTickLabel', {})
			end
            
            % Set remaining axes properties
% 			set(ax, ...
% 				'TickLength', [0.04 0.05], ...
% 				'TickDir', 'out', ...
% 				'LineWidth', 1.5, ...
% 				'Box', 'off', ...
% 				'FontSize', 8)
			set(ax, ...
				'TickLength', [0.02, 0.025], ...
				'TickDir', 'out', ...
				'LineWidth', 1, ...
				'box', 'off')
			
			% Make sure ax is current before returning
% 			axes(ax);
			
%             hha3=axes;
%             set(hha3,'color','none', ...
%                 'XTick',[],'XTickLabel',{},...
%                 'YTick',[],'YTickLabel',{},...
%                 'xlim', axTransformed,...
        end

        
		function [nelements , centers] = biexhist(Y, M, showPlot)
			%BIEXHIST Logiclly-scaled (Parks,et al.) histogram. 
			%   BIEXHIST(Y) bins the elements of Y into 10 equally spaced containers
			%   and produces a histogram plot of the results.  If Y is a
			%   matrix, hist works down the columns.
			% 
			%   N = hist(Y, M, showPlot), where M is a scalar, uses M bins.
			%	
			%	showPlot is an optional logical flag to plot the histogram.
			%
			%   Example:
			%       [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs('sample.fcs');
			%       GreenChannel = getChannel(fcshdr,'FIT');
			%       GreenData = fcsdat(:,GreenChannel);
			%       biexhist(greenData)
			%
			%   Written by
			%   Breanna Stillo
			%   bstillo@mit.edu
			%   Last Updated: 2014-10-31;
			%
			%   UPDATES:
			%   10/31 -- changed ylimits to reflect FlowJo.  Change default number of
			%   bins

			trf = Transforms();
			ylog = trf.lin2logicle(Y);
			
			if ~exist('M','var')
				M = round(77 * range(ylog));
			end
			[nelements, centers] = hist(ylog, M);
			
%             % Normalize to max
%             maxN = round(max(nelements) * 1e-3) / 1e-3;
% 			  if maxN < max(nelements)
% 				  maxN = maxN + 1000;
%             end            
%             nelements=nelements./max(nelements).*100;

			if (exist('showPlot', 'var') && showPlot)
				area(centers, nelements, 'FaceColor', [0.5 0.5 0.5], 'LineWidth', 1.5);
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
            
            % Check existence of optional input dataType, assign default
            %   faceColor and axPosition are also optional, but they are looked for later
            if (~exist('dataType', 'var'))
                dataType = 'raw';
            end
            
            % Check inputs
            checkInputs_histFit(data, channel, dataType);
            
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
            
            
            function checkInputs_histFit(data, channel, dataType)
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
            % Written by Ross Jones
            %   Weiss Lab, MIT
            %   2016-03-26
            %   
            % Update log:
            
            % Check inputs
            checkInputs_scatterBarPlot(data, channel, dataType);
            
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
            
            
            function checkInputs_scatterBarPlot(data, channel, dataType)
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
            % Written by Ross Jones
            %   Weiss Lab, MIT
            %   2016-03-26
            %   
            % Update log:
            
            % Check inputs
            checkInputs_batchBoxPlot(data, channel, dataType);
            
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
            
            
            function checkInputs_batchBoxPlot(data, channel, dataType)
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
            % Written by Ross Jones
            %   Weiss Lab, MIT
            %   2016-04-04
            %   
            % Update log:
            
            % Check inputs
            checkInputs_batchViolinPlot(data, channel, dataType);
            
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
            for i = 1:numel(data)
                
                % Assign x- and y-axis datasets
                ydata = Transforms.lin2logicle(data(i).(channel).(dataType));
                
                % Plot points
                Plotting.violinplot(ax, ydata, i, 'hist', colors(i, :));
            end
            
            % Set Y-axis to biexp
            Plotting.biexpAxes(ax, false, true);
            ax.XTick = 1:numel(data);
            ax.XTickLabelRotation = -45;
            ax.FontSize = 14;
            
            
            % --- Helper Function --- %
            
            
            function checkInputs_batchViolinPlot(data, channel, dataType)
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
            % Written by Ross Jones
            %   Weiss Lab, MIT
            %   2016-05-09
            %   
            % Update log:
            
            % Check inputs
            [data, channel, dataType] = checkInputs_lineDensityPlot(inputs);
            
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
						xdata = Transforms.lin2logicleMEF(xdata);
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
					Plotting.biexpAxesMEF(ax, true, false);
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
            
            
            function [data, channel, dataType] = checkInputs_lineDensityPlot(inputs)
                
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
        
		
		function singleLineDensity(ax, xdata, color, options)
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
			%						'counts'	Flag to plot bin counts rather than PDF
			%						'smooth'	Flag to smooth the histogram data
			
			checkInputs_singleLineDensity();
			
			% Fix negative infinite values by setting the resulting values to the 
			% otherwise minimum value.
			% Fix NaN and inf values by removing them.
			if any(isnan(xdata) | isinf(xdata))
				warning('Inf/NaN values detected - removing')
				valid = ~(isnan(xdata) | isinf(xdata));
				xdata = xdata(valid);
			end
			
			% Automatically find number of bins based on data
			if ismember('counts', options)
				normMode = 'count';
			else
				normMode = 'pdf';
			end
			[binCounts, edges] = histcounts(xdata, 30, 'Normalization', normMode);
			
			% Estimate bin centers by averaging the edges
			binCenters = zeros(numel(edges) - 1, 1);
			for j = 2:length(edges)
				binCenters(j - 1) = mean(edges([j - 1, j]));
			end
			
			% Ignore bins w/o data, otherwise gaps in the line plot will show up
			hasPoints = (binCounts > 0);
			binCounts = binCounts(hasPoints);
			binCenters = binCenters(hasPoints);
			if ismember('smooth', options)
				binCounts = smooth(binCounts);
			end
			
			% Plot w/ or w/o shading
			if ismember('shade', options)
				area(ax, binCenters, binCounts, 'FaceColor', color, 'linewidth', 1, 'FaceAlpha', 0.4);
				plot(ax, binCenters, binCounts, 'color', color, 'linewidth', 3);
				line(ax, [binCenters(1), binCenters(1)], [1, binCounts(1)], 'color', color, 'linewidth', 3);
				line(ax, [binCenters(end), binCenters(end)], [1, binCounts(end)], 'color', colors(i, :), 'linewidth', 3);
			else
				plot(ax, binCenters, binCounts, 'color', color, 'linewidth', 5);
			end
			
			
			% --- Helper Functions --- %
			
			
			function checkInputs_singleLineDensity()
				
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
				if ~exist('options', 'var')
					options = {};
				end
				validateattributes(options, {'cell', 'char'}, {}, mfilename, 'options', 4);
				
				% For simplicity, make options be a cell array
				if ischar(options), options = {options}; end
			end
		end
		
        
		function orthogonalityMatrix(A,labels)
			%orthogonalityMatrix(A,labels)
			%generates an orthogonality plot for the matrix A
			%   A = matrix of data vlaues (must be square)
			%   labels = labels for the matrix
			%
			%Ex:
			%   A =[ 4.5000    0.5000    0.2500    0.5500
			%     0.0500    3.5000    1.0000    0.2500
			%     0.5000    0.7500    4.2500    0.1500
			%     0.6500    0.4500    1.2500    4.7500];
			%   orthoMatrix(A,{'A1','A2','A3','A4'})
			%
			%   Written by
			%   Breanna Stillo
			%   bstillo@mit.edu
			%   Last Updated: 2014-10-14;


			map=[linspace(0,1)' zeros(100,1) linspace(1,0)'];
			colormap(map)
			imagesc(A)
			colorbar('southoutside')

			set(gca,'Xtick',1:length(A),...
				'Ytick',1:length(A),...
				'XAxisLocation','top')

            if exist('labels','var')
                set(gca,'XTickLabel',labels,'YTickLabel',labels)
            end
        end
		
        
        function [ax, h] = standardHeatmap(data, cmap, rowLabels, colLabels, norm, doCluster, symmetric, range, rowpdist, colpdist)
            % Generates a standard heatmap with fontsize 14 and data arranged nicely
            %
            %   Inputs:
            %
            %       data            The data to be plotted (matrix). Normally, heatmaps flip
            %                       the rows for some reason (like an image) - this is 
            %                       corrected such that the heatmap "starts" in the top left
            %       
            %       cmap            The colormap to use for the heatmap (Nx3 matrix of RGB colors)
            %
            %       rowLabels       Cell list of strings with row labels
            %
            %       colLabels       Cell list of strings with row labels (the function accounts 
            %                       for flipping as mentioned above). 
            %
            %       norm            (optional) The dimension to normalize on
            %                           1 or 'column' = within columns
            %                           2 or 'row' = within rows
            %                           3 or 'none' = no normalization
            %
            %       doCluster       (optional) Boolean value which flags to create a clustergram
            %                       rather than a heatmap. Clustergram is created with standard
            %                       euclidean distance calculation for both rows and cols.
            %
            %       symmetric       (optional) Boolean value which flags the heatmap to be 
            %                       symmetric around zero. (Default = False)
            %   Outputs: 
            %
            %       ax              A handle to the figure axes created
            %       h               A handle to the heatmap/clustergram object
            % 
            % Written by Ross Jones
            % jonesr18@mit.edu
            % Weiss Lab, MIT
            % 2016-04-13
            %
            % Update Log:
            %   2016-05-02 Added symmetric argument
            
            % Invert data rows because stupid heatmap
            flippedData = flipud(data);
            
            % Create heatmap/clustergram
            if (exist('doCluster', 'var') && doCluster)
                h = clustergram(flippedData);
                if exist('rowpdist', 'var')
                    h.RowPDist = rowpdist;
                else
                    h.RowPDist = 'correlation';
                end
                if exist('colpdist', 'var')
                    h.ColumnPDist = colpdist;
                else
                    h.ColumnPDist = 'correlation';
                end
            else
                doCluster = false;
                h = HeatMap(flippedData);
            end
            if exist('norm', 'var')
                h.Standardize = upper(norm);
                if exist('range', 'var')
                    h.DisplayRange = range; 
                end
            end
            if exist('symmetric', 'var')
                h.Symmetric = symmetric;
            else
                h.Symmetric = false;
            end
            h.Colormap = cmap;
            h.RowLabels = fliplr(rowLabels);
            h.ColumnLabels = colLabels;
            
            h.view();
            
            % Create figure/axes object and update font size / add colorbar
            ax = h.plot();
            if (~doCluster)
                % The colorbar throws off the dendrogram, so it needs to be generated in other ways
                colorbar('peer', ax, 'EastOutside');
            end
            ax.FontSize = 14;
        end
	end
	
end