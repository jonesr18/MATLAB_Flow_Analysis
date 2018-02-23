classdef Compensation < handle
	% Compilation of compensation methods from the Weiss Lab flow cytometry MATLAB repository
	%
	%	Methods are implemented as static functions, so they can be called directly
	%	or after creating a handle to this class.
	%
	%		Ex:
	%		data = Compensation.openFiles(TF_marker, args);
	%
	%		Ex2:
	%		comp = Compensation();
	%		compData = comp.compensateMatrix(inputData, coefficients);
	%
    %   Functions:
    %
    %       [xBleed, yComp] = compensateSC(JCFile, singleColorFile, correctFile, bC, cC, plotsOn)
    %       compData		= compensateBatchSC(jcData, scData, scChannels, inputData, fixChannels, dataType, gate, plotsOn)
    %		[data, wtData, scData, fitParams] = compensateMatrixBatch(data, channels, wtData, scData, dataType, gate, plotsOn)
	%       compdData		= compensateMatrix(inputData, coefficients)
    %       YI				= lsq_lut_piecewise(x, y, XI)
    %
    % Written/Compiled by Ross Jones
    % Weiss Lab, MIT
    % Last updated 2017-09-26
	
	methods (Static)
		
		function [xBleed, yComp] = compensateSC(JCFile, singleColorFile, correctFile, bC, cC, plotsOn)
			% [xBleed,yComp] = compensateSC(JCFile, singleColorFile, correctFile, bC, cC,plotsOn)
			% single-color compensation
			%   JCFile = Just Cell File
			%   singleColorFile = FCS file for single color control
			%   correctFile = FCS file that needs to be corrected
			%   bC = channel number of the color that is bleeding and needs to be
			%        corrected for in other channels
			%   cC = channel number of the color that needs to be corrected
			%   plotsOn = optional. Turns plots on for visualization
			%     
			%   Example:
			%       [fcsdat, ~, ~, ~] = fca_readfcs('sample.fcs');
			%       GreenChannel = getChannel(fcshdr,'FIT');
			%       RedChannel = getChannel(fcshdr,'Red');
			%       GreenData = fcsdat(:,GreenChannel);
			%       RedData = fcsdat(:,RedChannel);
			%       figure()
			%       biexplot(redData,greenData)
			%       [newRed,newGreen]=compensateSC('JCsample.fcs','RedSample.fcs','sample.fcs',RedChannel,GreenChannel);
			%       figure()
			%       biexplot(newred,newgreen)
			%
			%   Written by
			%   Breanna Stillo
			%   bstillo@mit.edu
			%   Last Updated: 2014-10-21;

			[fcsdat, ~, ~, ~] = fca_readfcs(JCFile);
			cValsJC = fcsdat(:, cC);
			median_cJC = median(cValsJC);
			
			[fcsdat, ~, ~, ~] = fca_readfcs(singleColorFile);
			bValsSC = fcsdat(:, bC);
			cValsSC = fcsdat(:, cC);

			XI = [0 10^2 10^2.5 10^3 10^3.5 10^4 10^4.5 max(bValsSC)];
			YI = Compensation.lsq_lut_piecewise(bValsSC, cValsSC, XI);
			
			if exist('plotsOn','var')
				if plotsOn
					figure()
					biexplot(bValsSC, cValsSC)
					hold on
					biexplot(XI, YI, 'r-')
				end
			end
			
			[fcsdat, ~, ~, ~] = fca_readfcs(correctFile);
			bValsCF = fcsdat(:, bC);
			cValsCF = fcsdat(:, cC);
			
			yerror = compEval(bValsCF);
			cReal = cValsCF - yerror + median_cJC;
			
			xBleed = bValsCF;
			yComp = cReal;
			
			
			% --- Helper Function --- %
			
			
			function y = compEval(x)
						
				yBuff = zeros(length(x), 1);
				for w = 1:length(x)
					eqn = find(XI > x(w), 1, 'first');
					if isempty(eqn)
						eqn = length(XI);
					end
					if eqn == 1
						eqn = 2;
					end

					a = (YI(eqn) - YI(eqn - 1)) ./ (XI(eqn) - XI(eqn - 1));
					b = YI(eqn) - a * XI(eqn);
					yBuff(w) = a * x(w) + b;
				end
				
				y = yBuff;
					
			end
			
        end
        
        
        function compData = compensateBatchSC(jcData, scData, scChannels, inputData, fixChannels, dataType, gate, plotsOn)
            % compensateBatchSC compensates data using a fitting routine to single-color bleeding.
            % 
            %   compData = compensateBatchSC(jcData, scData, scChannels, inputData, fixChannels, gate, plotsOn)
            %   
            %       Inputs: 
            %           jcData          A standard struct with untransfected cells data
            %            (struct)       
            %
            %           scData          A standard struct with single color samples data
            %            (struct)       
            %
            %           scChannels      A cell array of single color channels - entries should 
            %            (cell)         be strings with the names given by Fortessa
            %                               eg: {'PE_Texas_Red_A', 'FITC-A'}
            %                               Note this is case-sensitive
            %                           The channels must match with the samples in scData!
            %                               
            %             ****              The order MUST BE THE SAME as the single-color samples
            %                               given in scData. This will not be auto-checked!
            %
            %           inputData       A stardard struct with channel data to be compensated. 
            %            (struct)       'fixChannels' should all exist in this struct.
            %
            %           fixChannels     A cell array of strings indicating the channels to be fixed
            %            (cell)         by the compensation procedure.
            %
			%			dataType (opt)	A string indicating which data type to compensate with. 
			%			 (string)		
            %
            %           gate (opt)      A string indicating what gated population to compensate
            %            (string)       with, for example, 'P3'. Must match entries in 'gates'
            %                           field of all data structs.
            %
            %           plotsOn (opt)   A boolean indicating whether to plot or not
            %            (logical)
            %
            %       Outputs:
            %           compData        A copy of the inputData struct with added scComp field
            %            (struct)       for each compensated channel.
            %
            %   Modified from Bre's compensateSC method above for automated batch processing and
            %   a simpler interface for the user using my data structure.
            %
            %   Ross Jones
            %   Weiss Lab
            %   2015-11-17
            
            % Check the inputs
            zCheckInputs_compensateBatchSC();
                        
            % Allocate output structure
            fields = fieldnames(inputData);
            nSamplesToComp = length(inputData);
            compData(nSamplesToComp) = cell2struct(cell(length(fields), 1), fields, 1);
            
            % Compensate using each channel pair
            for i = 1:nSamplesToComp
                compDataSample = inputData(i);
                for j = 1:numel(scChannels)
                    bleedChannel = scChannels{j};
                    for k = 1:numel(fixChannels)

                        % Skip compensation for same channel matches
                        if strcmpi(scChannels{j}, fixChannels{k})
                            continue
                        end
                        
                        % Process compensation and store compensated data in output struct
                        % Note that on each call of proscesCompensation, the variable 
                        % compDataSample is updated with compensation data. Thus, this data
                        % is added to the function output after each channel has been processed
                        
                        fixChannel = fixChannels{k};
                        
                        % -- It turns out this is slower to load/save the YI vector - computation is
                        %    relatively fast!
                        processCompensation(bleedChannel, fixChannel, plotsOn); 

                    end
                end
                compData(i) = compDataSample;
            end
            
            
            % --- Helper functions --- %
            
            
            function zCheckInputs_compensateBatchSC()
                
                % Check types
                validateattributes(jcData, {'struct'}, {}, mfilename, 'jcData', 1);
                validateattributes(scData, {'struct'}, {}, mfilename, 'scData', 2);
                validateattributes(scChannels, {'cell', 'char'}, {}, mfilename, 'channels', 3);
                validateattributes(inputData, {'struct'}, {}, mfilename, 'inputData', 4);
                validateattributes(fixChannels, {'cell', 'char'}, {}, mfilename, 'channels', 5);
				if exist('dataType', 'var')
					validatestring(dataType, fieldnames(jcData(1).(scChannels{1})), mfilename, 'dataType', 5)
				end
                if exist('gate', 'var')
                    validatestring(gate, fieldnames(jcData(1).gates), mfilename, 'gate', 6);
                end
                if exist('plotsOn', 'var')
                    validateattributes(plotsOn, {'logical'}, {}, mfilename, 'plotsOn', 7);
                else
                    plotsOn = false;
                end
                
                % Check if number of channels is less than or equal to length of single-color struct
                %   (should be 1:1 matching, or fewer if not all are needed)
                if ischar(scChannels), scChannels = {scChannels}; end
                assert(length(scChannels) <= length(scData), ...
                    'Number of channels (%d) is greater than the number of single-color controls (%d)', ...
                    length(scChannels), length(scData));
                
                % Check if each given single color channel is in the measured channels
                for iii = 1:numel(scChannels)
                    channel = scChannels{iii};
                    if isempty(intersect(channel, fieldnames(scData(iii))))
                        error('Given scChannel: %s is not valid', channel{:});
                    end
                end
                                
                % Check if each given fix channel is in the measured channels
                if ischar(fixChannels), fixChannels = {fixChannels}; end
                for channel = reshape(fixChannels, 1, [])
                    if isempty(intersect(channel{:}, fieldnames(inputData(1))))
                        error('Given fixChannel: %s is not valid', channel{:});
                    end
                end
            end
            
            
            function YI = processCompensation(bleedChannel, fixChannel, plotsOn)
                % Avoid overhead of repeatedly passing large data structures as inputs or 
                % reloading data by using global variables. 
                
                
                % Get gate indexes if applicable
                if exist('gate', 'var')
                    jcGate = jcData.gates.(gate);
                    scGate = scData(j).gates.(gate);
                    inGate = compDataSample.gates.(gate);
                end
                
                % Extract median value of untransfected cells for the fix channel
                jcMedianFixChannel = median(jcData.(fixChannel).raw(jcGate));
                
                % Extract single color data for bleed and fix channels
                %   Sub index with j, which is the scFile corresponding with bleedchannel
                scBleedChannel = scData(j).(bleedChannel).raw(scGate);
                scFixChannel = scData(j).(fixChannel).raw(scGate);
                
                % Extract input data for the bleed and fix channels
                cdsBleed = compDataSample.(bleedChannel);
                if isfield(cdsBleed, 'scComp')
                    inBleedChannel = cdsBleed.scComp;
                else
                    inBleedChannel = cdsBleed.raw(inGate);
                end
                cdsFix = compDataSample.(fixChannel);
                if isfield(cdsFix, 'scComp')
                    inFixChannel = cdsFix.scComp;
                else
                    inFixChannel = cdsFix.raw(inGate);
                end
                
                % Setup arrays for curve fitting
                XI = [-10^4, -10^3, -10^2 0, 10^2, 10^2.5, 10^3, 10^3.5, 10^4, 10^4.5, 10^5];
                YI = Compensation.lsq_lut_piecewise(scBleedChannel, scFixChannel, XI);
                
                % Evaluate bleeding: subtract the difference between the yerror (the amount of
                % bleed as determined by the single color controls) and the background.
                yerror = interp1(XI, YI, inBleedChannel);
                try
                    compdInFixChannel = inFixChannel - (yerror - jcMedianFixChannel);
                catch ME
                    whos inFixChannel yerror 
                    cdsBleed %#ok<NOPRT>
                    cdsFix %#ok<NOPRT>
                    error(ME)
                end
                compDataSample.(fixChannel).scComp = compdInFixChannel;
%                 valid = ~(isnan(compdInFixChannel) | isinf(compdInFixChannel));
%                 compDataSample.(fixChannel).scComp = compdInFixChannel(valid);
%                 if isfield(cdsBleed, 'scComp')
%                     compDataSample.(bleedChannel).scComp = inBleedChannel(valid);
%                 end
%                 compDataSample.gates.(gate) = valid;
                
                % Test plotting
                if exist('plotsOn','var') && plotsOn
                    figure()
                    biexplot(scBleedChannel, scFixChannel, {'.', 'markersize', 5})
                    hold on
                    biexplot(XI, YI, '-')
                    hold on
                    biexplot(inBleedChannel, yerror, {'.', 'markersize', 5})
                    hold on
                    biexplot(inBleedChannel, inFixChannel, {'.', 'markersize', 5})
                    hold on
                    biexplot(inBleedChannel, compdInFixChannel, {'.', 'markersize', 5})
                    pause()
                end

            end
		end
		
		
		function [coeffs, ints, figFits] = computeCoeffs(scData, channels, options)
			% Finds the coefficients for linear fits between each channel
			% when observing bleed-through.
			%
			%	[coeffs, ints, figFits] = computeCoeffs(dataMatrix, channels, plotsOn, minFunc)
			%
			%	Inputs
			%
			%		scData		<cell> A cell list of C NxC data matrices containing 
			%					N cell fluroescenct values in C channels from C 
			%					single-color controls. The ordering of the columns 
			%					(channels) should match the ordering of the 
			%					corresponding controls in the cell list. 
			%
			%		channels	<cell> A cell list of C channel names
			%					corresponding with the cell list and
			%					matrix data in scData.
			%
			%		options		<struct> (optional) Optional property-value pairs:
			%						'minFunc':	A user-defined function for
			%									residual minimzation during
			%									fitting. Default = @(x) x.^2
			%									(least-squares approximation)
			%						'plotsOn':	If TRUE, shows the compensation
			%									plots, which are passed back to
			%									the caller (otherwise returns empty).
			%						'doMEF':	If TRUE, does logicle conversion
			%									with MEF-unit scaling
			%						'params':	<params> enables setting the
			%									logicle function parameters
			%									(see Transforms.lin2logicle())
			%
			%	Outputs
			%
			%		coeffs		<numeric> A CxC matrix of linear coefficients
			%					representing the bleed-through between channels
			%
			%		ints		<numeric> A Cx1 vector of linear intercepts
			%					representing autofluorescence in each channel
			%
			%		figFits		<struct> If options.plotsOn = TRUE, this returns 
			%					a handle to the generated figure with overlaid
			%					fit lines on the fitted data. If options.plotsOn
			%					is not passed or is FALSE, this returns empty.
			%
			%	Implementation notes
			%
			%		The function simultaneously fits coeffs and ints to all the
			%		data from each control and channel. ints is fixed for each
			%		channel, so some fits will appear off before compensation is
			%		applied, after which they are corrected. 
			%
			%		We also discard any value in a bleed channel that is > 10x
			%		the average value, as these often are spillover between
			%		tubes during data collection and can mess up compensation. 
			%
			% Update log:
			% TODO : Add back the outlier exclusion
			
			% Check inputs
			zCheckInputs_computeCoeffs();
			
			% Exclude outliers by ignoring any point that is >10x greater 
			% on the fix axis as the bleed axis. 
			for chB = 1:numel(channels)
				nonChB = setdiff(1:numel(channels), chB);
				outliers = any(scData{chB}(:, nonChB) > (10 * abs(scData{chB}(:, chB))), 2);
				scData{chB} = scData{chB}(~outliers, :);
			end
			
			% Equalize the number of points between each control
			numPoints = min(cellfun(@(x) size(x, 1), scData));
			scData = cellfun(@(x) x(FlowAnalysis.subSample(size(x, 1), numPoints), :), ...
				scData, 'uniformoutput', false);
			
			% Fitting initial conditions
			A0 = 10 * ones(numel(channels), 1);					% Intercepts (Autofluorescence)
			K0 = zeros(numel(channels)^2 - numel(channels), 1);	% Coefficients (Bleed-through)
			
			% Iterate over each pair of channels and compute fits
			optimOptions = optimset('Display', 'off');
			[minResult, fval] = fminsearch(@(x) fitFunc(x, options.minFunc), [A0; K0], optimOptions);
			fprintf('Linear fits obtained with obj func val %.2f\n', fval);
			ints = minResult(1:numel(channels));
			coeffs = eye(numel(channels));
			coeffs(~logical(eye(numel(channels)))) = minResult(numel(channels) + 1 : end);
			
			% Set up figure to view fitting (if applicable)
			if ~options.plotsOn
				figFits = [];
				return % No need to process the rest of the code
			end
			
			% Prep figure
			figFits = figure();
			spIdx = 0;
			xrange = logspace(0, ceil(log10(max(cellfun(@(x) max(x(:)), scData)))), 100);
			
			for chF = 1:numel(channels)
				for chB = 1:numel(channels) 
					
					fitVals = xrange * coeffs(chF, chB) + ints(chF);
					
					spIdx = spIdx + 1;
					ax = subplot(numel(channels), numel(channels), spIdx);
					hold(ax, 'on')
					
					plot(ax, Transforms.lin2logicle(scData{chB}(:, chB), ...
								 options.doMEF, options.logicle), ...
							 Transforms.lin2logicle(scData{chB}(:, chF), ...
								 options.doMEF, options.logicle), ...
							 '.', 'MarkerSize', 4)
					plot(ax, Transforms.lin2logicle(xrange, ...
								 options.doMEF, options.logicle), ...
							 Transforms.lin2logicle(fitVals, ...
								 options.doMEF, options.logicle), ...
							 '-', 'linewidth', 4)
					Plotting.biexpAxes(ax, true, true, false, ...
							 options.doMEF, options.logicle);
					
					% Axis labeling
					title(sprintf('Slope: %.2f | Intercept: %.2f', ...
						coeffs(chF, chB), ints(chF)), 'fontsize', 14)
					if (chF == numel(channels))
						xlabel(strrep(channels{chB}, '_', '-'))
					end
					if (chB == 1)
						ylabel(strrep(channels{chF}, '_', '-'))
					end
				end
			end
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_computeCoeffs()
				
				validateattributes(scData, {'cell'}, {}, mfilename, 'scData', 1);
				for sc = 1:numel(scData) 
					validateattributes(scData{sc}, {'numeric'}, {}, mfilename, sprintf('scData{%d}', sc), 1);
				end
				validateattributes(channels, {'cell'}, {}, mfilename, 'channels', 2);
				assert(numel(channels) == numel(scData), 'Incorrect number of channel controls or labels!')
				assert(numel(channels) == size(scData{1}, 2), 'Incorrect number of channel data or labels!');
				assert(numel(channels) > 1, 'Compensation with just one channel is useless!');
				
				if ~exist('options', 'var'), options = struct(); end
				if ~isfield(options, 'plotsOn'), options.plotsOn = false; end
				options.plotsOn = all(logical(options.plotsOn));
				if isfield(options, 'minFunc')
					validateattributes(options.minFunc, {'function_handle'}, {}, mfilename, 'options.minFunc', 4)
				else
					options.minFunc = @(x) sum(x.^2);
				end
				if ~isfield(options, 'logicle'), options.logicle = struct(); end
				if ~isfield(options, 'doMEF'), options.doMEF = false; end
				
			end
			
			
			function out = fitFunc(p, minFunc)
				% Fit function for computing linear fits between bleed 
				% and fix channels in each scData control
				
				% Think about trying to fit everything at once?
				
				% Setup fit matrix
				A = p(1:numel(channels));
				K = eye(numel(channels));
				K(~logical(eye(numel(channels)))) = p(numel(channels) + 1 : end);
				
				% Find regression with each channel
				residuals = zeros(1, numel(channels));
				for ch = 1:numel(channels)
					fixedData = K \ (scData{ch}' - A);
					chans = setdiff(1:numel(channels), ch); % Skip ch since is not min to 0
					residuals(ch) = minFunc(reshape(fixedData(chans, :), 1, []));
				end
				
				out = sum(residuals);
			end
		end
		
		
        function compData = matrixComp(uncompData, coeffs)
            % Compensates data using linear coefficients for bleed-through between channels.
            % 
            %   compData = matrixComp(uncompData, coeffs)
            %   
            %   Inputs: 
            %       uncompData		A CxN matrix of data points corresponding with
			%						N cells in C channels to be compensated. 
            %
            %       coeffs			A CxC matrix of coefficients corresponding with
            %						the slope of lines of best fit between C channels 
            %						bleeding into each other. 
            %
            %       Outputs:
            %           compData	An NxM matrix with compensated data.
            %   
            %   This method solves the linear set of equations:
            %       X_observed = A * X_real
            %
            %   'uncompData' corresponds with 'X_observed' and 'coeffs'
            %   corresponds with the matrix 'A'. We invert 'A' to solve
            %   for X_real, which is returned as 'compdData'.
            %
            % Ross Jones
            % Weiss Lab
            % 2017-06-06
            %
            % Update Log:
            %
            % ...
            
            % Check inputs are valid
            zCheckInputs_compensateMatrix();
            
            % Compensate data
            compData = coeffs \ uncompData;
            
            
            % --- Helper functions --- %
            
            
            function zCheckInputs_compensateMatrix()
                
                % Check types
                validateattributes(uncompData, {'numeric'}, {}, mfilename, 'uncompData', 1);
                validateattributes(coeffs, {'numeric'}, {}, mfilename, 'coeffs', 2);
                
                % Check sizes
                assert(size(uncompData, 1) == size(coeffs, 2), ...
                    '# Channels in ''uncompData'' (%d) and ''coeffs'' (%d) are not the same!', ...
                    size(uncompData, 1), size(coeffs, 2))
                assert(size(coeffs, 1) == size(coeffs, 2), ...
                    '''coeffs'' is not square! (H = %d, W = %d)', ...
                    size(coeffs, 1), size(coeffs, 2))
            end
        end
		
		
		function YI = lsq_lut_piecewise(x, y, XI)
			% LSQ_LUT_PIECEWISE Piecewise linear interpolation for 1-D interpolation (table lookup)
			%   YI = lsq_lut_piecewise( x, y, XI ) obtain optimal (least-square sense)
			%   vector to be used with linear interpolation routine.
			%   The target is finding Y given X the minimization of function 
			%           f = |y-interp1(XI,YI,x)|^2
			%   
			%   INPUT
			%       x measured data vector
			%       y measured data vector
			%       XI break points of 1-D table
			%
			%   OUTPUT
			%       YI interpolation points of 1-D table
			%           y = interp1(XI,YI,x)
			%
            
            % Turn off this warning, which results from using very negative values in XI
            warning('off', 'MATLAB:rankDeficientMatrix')
            
            % Check vector sizes
			if size(x, 2) ~= 1
				error('Vector x must have dimension n x 1.');   
			elseif size(y, 2) ~= 1
				error('Vector y must have dimension n x 1.');    
			elseif size(x, 1) ~= size(y, 1)
				error('Vector x and y must have dimension n x 1.'); 
			end

			% matrix defined by x measurements
			A = sparse([]); 

			% vector for y measurements
			Y = []; 

			for j = 2:length(XI)
				
				% get index of points in bin [XI(j-1) XI(j)]
				ix = (x >= XI(j - 1) & x < XI(j) );
				
				% check if we have data points in bin
                % doesn't matter
% 				if ~any(ix)
% 					warning('Bin [%f %f] has no data points, check estimation. Please re-define X vector accordingly.',XI(j-1),XI(j));
% 				end
				
				% get x and y data subset
				x_ = x(ix);
				y_ = y(ix);
				
				% create temporary matrix to be added to A
				tmp = [(1 - (x_ - XI(j - 1)) / (XI(j) - XI(j - 1))),...
                       (    (x_ - XI(j - 1)) / (XI(j) - XI(j - 1))) ];
				
				% build matrix of measurement with constraints
				[m1, n1] = size(A);
				[m2, n2] = size(tmp);
				A = [A, zeros(m1, n2 - 1); 
                     zeros(m2, n1 - 1), tmp]; %#ok<AGROW>
				
				% concatenate y measurements of bin
				Y = [Y; y_]; %#ok<AGROW>
			end

			% obtain least-squares Y estimation
			YI = (A \ Y)';

		end
		
	end
	
end