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
				switch mode %#ok<ALIGN>
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
						nColors = 1000;
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
						nColors = 1000;
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
				
                nColors = 1000;
                cm = colorMap.getColormap(nColors);
				MIN = min(sortedZ);
					if (nonZero), MIN = max(0, MIN); end
                colors = cm(ceil((sortedZ - MIN) ./ (4.5 - MIN) .* (nColors - 1) + 1), :);
            end
            
            % Plot data in sorted order so highest density points are plotted last.
            scatter(ax, x(sortIdx), y(sortIdx), 10, colors, 'filled')
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
            checkInputs(ydata, xcenter, mode, faceColor); 
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
            
            
            function checkInputs(ydata, xcenter, mode, faceColor)
              
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
		
        
        function biexpAxes(ax, biexpX, biexpY, textOn)
            % Tranforms the given axes to a scale used for biexponential veiws.
            %   
            %   Inputs
            %       
            %       biexpX (logical)        TRUE - make x-axis biexponential
            %       biexpY (logical)        TRUE - make y-axis biexponential
            %       textOn (logical)        TRUE - include axis text
            %
            %   If inputs are not given, default to biexponential for both X 
            %   and Y axes with axis labels on.
            
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
            
            % Set values
			AXES_MAX = 2^18;
			AXES_MIN = -150;
            axTransformed = Transforms.lin2logicle([AXES_MIN, AXES_MAX]);

            % Set axes properties
            largeTickVals = Transforms.lin2logicle([-1e2 -1e1 -1e0 0 1e0 1e1 1e2 1e3 1e4 1e5 ]);
            set(ax, ...
               'ticklength', [0.04 0.05], ...
               'TickDir', 'out', ...
               'LineWidth', 1.5, ...
               'box', 'off', ...
               'fontsize', 8)
            
            % Write axes labels
            if (exist('textOn', 'var') && ~textOn)
                axTxt = '';
            else
                axTxt = {'-10^2', '', '', '  0^{ }', '', '', ' 10^2', ' 10^3', ' 10^4', ' 10^5'};
            end
                        
            if (biexpX)
                set(ax, ...
                   'xlim', axTransformed, ...
                   'xtick', largeTickVals, ...
                   'xticklabel', axTxt) 
            end
            
            if (biexpY)
                set(ax, ...
                   'ylim', axTransformed, ...
                   'ytick', largeTickVals, ...
                   'yticklabel', axTxt)
            end
            
%             % Add smaller ticks
%             tickVals = sort( ...
%                 [-10^2, -(1:9).*10^1, -(1:9).*10^0, 0, ...
%                  (1:9).*10^0, (1:9).*10^1, (1:9).*10^2, (1:9).*10^3, ...
%                  (1:9).*10^4, 1e5]);
%             set(ax, ...
%                 'color', 'none', ...
%                 'XTick', lin2logicle(tickVals), ...
%                 'XTickLabel', {}, ...
%                 'YTick', lin2logicle(tickVals), ...
%                 'YTickLabel', {}, ...
%                 'xlim', lin2logicle([AXES_MIN, AXES_MAX]), ...
%                 'TickDir', 'out', ...
%                 'ylim', lin2logicle([AXES_MIN, AXES_MAX]), ...
%                 'LineWidth', 1)
% 
%             hha3=axes;
%             set(hha3,'color','none', ...
%                 'XTick',[],'XTickLabel',{},...
%                 'YTick',[],'YTickLabel',{},...
%                 'xlim',lin2logicle([AXES_MIN, AXES_MAX]),...
%                 'ylim',lin2logicle([AXES_MIN, AXES_MAX]),'box','on','LineWidth',1.5) 

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
            checkInputs(data, channel, dataType);
            
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
            
            
            function checkInputs(data, channel, dataType)
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
            checkInputs(data, channel, dataType);
            
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
            
            
            function checkInputs(data, channel, dataType)
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
            checkInputs(data, channel, dataType);
            
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
            
            
            function checkInputs(data, channel, dataType)
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
            checkInputs(data, channel, dataType);
            
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
            
            
            function checkInputs(data, channel, dataType)
                % Checks the inputs to make sure they are valid
                validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                validatestring(channel, fieldnames(data(1)), mfilename, 'channel', 2);
                validatestring(dataType, fieldnames(data(1).(channel)), mfilename, 'dataType', 3);
            end
        end
        
        
        function lineDensityPlot(inputs)
            % Creates a line density plot for each element in the given data struct. Lines are
            % essentially the connected quantities of bins from a histogram of a channel's data. 
            %
            %   Inputs (struct array w/ fields below):
            %       ax              (optional) The axes to plot on
            %       data            A standard struct with the data to be plotted
            %                        - Must have 'channel' and 'dataType' given by user
            %                        - elements are plotted in a row in order
            %       channel         The channel to plot data from
            %       dataType        The data type to use, eg 'raw', or 'scComp'
            %       colors          (optional) An Nx3 array of rbg color triplets - determines 
            %                       the color order used for the face of each violin
            %                        - N must be numel(data)
            %                        - Defaults to standard MATLAB sequence if no input given
            %       xscale          (optional) String: 'log', 'logicle'. Determines the scaling
            %                       for the x-axis. Default = logicle.
            %       shade           (optional) Logical: Determines whether to shade in the area
            %                       under the line plot. Default = false.
            %
            % Written by Ross Jones
            %   Weiss Lab, MIT
            %   2016-05-09
            %   
            % Update log:
            
            % Check inputs
            [data, channel, dataType] = checkInputs(inputs);
            
            % Ensure axis is held on
            if (~isfield(inputs, 'ax'))
                ax = gca();
            else
                ax = inputs.ax;
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
            
            colors
            
            % Plot data
            for i = 1:numel(data)
                
                % Assign x dataset
                if isfield(inputs, 'xscale') 
                    if strcmpi(inputs.xscale, 'linear')
                        xdata = data(i).(channel).(dataType);
                    elseif strcmpi(inputs.xscale, 'log')
                        xdata = log10(data(i).(channel).(dataType));
                    elseif strcmpi(inputs.xscale, 'logicle')
                        xdata = Transforms.lin2logicle(data(i).(channel).(dataType));
                    else
                        error('Unrecognized xscale given: %s', inputs.xscale);
                    end
                else
                    xdata = Transforms.lin2logicle(data(i).(channel).(dataType));
                end
                x = real(xdata);
                
                % Fix negative infinite values by setting the resulting values to the 
                % otherwise minimum value.
                % Fix NaN and inf values by removing them.
                if any(isnan(x) | isinf(x))
                    warning('Inf/NaN values detected - removing')
                    valid = ~(isnan(x) | isinf(x));
                    x = x(valid);
                end
                
                % Automatically find number of bins based on data
                [binCounts, edges] = histcounts(x, 30);
                
                % Estimate bin centers by averaging the edges
                binCenters = zeros(numel(edges) - 1, 1);
                for j = 2:length(edges)
                    binCenters(j - 1) = mean(edges([j - 1, j]));
                end
                
                % Ignore bins w/o data, otherwise gaps in the line plot will show up (Y log scale)
                hasPoints = (binCounts > 0);
                binCounts = binCounts(hasPoints);
                binCenters = binCenters(hasPoints);
                if (isfield(inputs, 'xscale') && strcmpi(inputs.xscale, 'log'))
                    binCenters = 10.^binCenters;
                end
                
                if (isfield(inputs, 'shade') && inputs.shade)
                    area(ax, binCenters, binCounts, 'FaceColor', colors(i, :), 'linewidth', 1, 'FaceAlpha', 0.4);
                    plot(ax, binCenters, binCounts, 'color', colors(i, :), 'linewidth', 3);
                    line(ax, [binCenters(1), binCenters(1)], [1, binCounts(1)], 'color', colors(i, :), 'linewidth', 3);
                    line(ax, [binCenters(end), binCenters(end)], [1, binCounts(end)], 'color', colors(i, :), 'linewidth', 3);
                else
                    plot(ax, binCenters, binCounts, 'color', colors(i, :), 'linewidth', 5);
                end
            end
            
            % Scale axes
            if isfield(inputs, 'xscale')
                if strcmpi(inputs.xscale, 'log')
                    ax.XScale = 'log';
                elseif strcmpi(inputs.xscale, 'linear')
                    ax.XScale = 'linear';
                end 
            else
                Plotting.biexpAxes(ax, true, false);
            end
            ax.YScale = 'log';
            ax.FontSize = 16;
             
            
            % --- Helper Function --- %
            
            
            function [data, channel, dataType] = checkInputs(inputs)
                
                % Ensure necessary inputs are present
                reqFields = {'data', 'channel', 'dataType'};
                hasFields = isfield(inputs, reqFields);
                if ~all(hasFields)
                    error('Missing field: %s\n', reqFields{~hasFields})
                end
                
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