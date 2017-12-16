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
            checkInputs_compensateBatchSC();
                        
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
            
            
            function checkInputs_compensateBatchSC()
                
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
        
        
        function [sampleData, controlData, fitParams] = ...
                        compensateMatrixBatch(sampleData, controlData, channels, dataType, gate, plotsOn)
            % compensateMatrixBatch compensates an entire standard data struct by generating 
            % coefficients with linear fitting on single-color control data and processing the
            % data to be sent piecewise to Compensation.compensateMatrix(). 
            % 
            % The method first subtracts autofluorescence based on the untransfected cell
            % sample provided (see below)
            %
            %   compData = compensateMatrixBatch(sampleData, controlData, channels, dataType, gate, plotsOn)
            %   
            %       Inputs: 
            %           sampleData      A stardard struct with channel data to be compensated. 
            %            (struct)       'channels' should all exist in this struct.
            %           
			%			controlData     A standard struct with control data to help with
			%			 (struct) 		compensation and to be compensated. 
            %							
            %                       *** The order of controls is very particular:
			%							The first N elements in the struct array 
			%							are single-color controls used for calculating
			%							bleed-through, and must be in the same order as 
			%							their corresponding colors in 'channels'
			%							-- This will not be auto-checked!
			%							
			%							The very last element must correspond with the 
			%							wild-type cell data (untransfected).
			%							
			%							Between the WT data and SC data, there
			%							can be two-color or other controls which
			%							are not used for compensation. 
			%
            %           channels        A cell array of strings indicating the channels to be fixed
            %            (cell)         by the compensation procedure.
            %
            %                           Entries should be strings with the names given by Fortessa
            %                           but with spaces and dashes replaced by underscores
            %                               eg: {'PE_Texas_Red_A', 'FITC_A'}
            %                               Note: this is case-sensitive
			%
			%			dataType (opt)	A string indicating which data type to compensate with. 
			%			 (string)		
            %
            %           gate (opt)      A string indicating which gated population to compensate
            %            (string)       with, for example, 'P3'. Must match entries in 'gates'
            %                           field of all data structs.
            %
            %           plotsOn (opt)   A boolean indicating whether to plot or not
            %            (logical)
            %
            %       Outputs:
            %           sampleData      A copy of the data struct with added 'afs' and 'mComp' 
            %            (struct)       fields for each compensated channel.
            %                               afs = autofluorescence subtracted
            %                               mComp = matrix compensated
            %
            %           controlData     A copy of the controlData struct with added 'afs' and 
            %            (struct)       'mComp' fields for each compensated channel.
			%
			%			fitParams		(Optional) The parameters fit during compensation
            %
            %
            % Ross Jones, Weiss Lab
            % jonesr18@mit.edu
            % Created 2017-06-06
            % 
            % Update Log: 
            %
            %	2017-07-16		Added dataType input
			%	2017-08-22		Added full compensation parameter output,
			%					changed autofluorescence calculation
			%	2017-10-03		Changed interface to accept single-struct controlData
			%
			%	...
			%
			% Autofluorescence:		[ 4,  2, 0,   1]
			% Summed Intercepts:	[13, 11, 0, -15]
			% Ints after AF sub:	[ 1,  5, 0, -18]
			% AF Sub effect:	  - [12,  6, 0,   3] = -3 * AF
            
            % Check the inputs - makes no changes to data structures
            checkInputs_compensateMatrixBatch();
            
            % Find medians for each channel in jcData
% 			autofluor = getAutofluor();
            
            % Find coefficients using least-squares linear fits - uses afs data
            [coeffs, ints] = getCoefficients();
			
			% Calculate autofluorescence from scData fit Y-intercepts
			autofluor = sum(ints, 2) / (numel(channels) - 1); 
% 			autofluor = zeros(numel(channels), 1);
			
			% Subtract autofluorescence - adds autofluorescence subtracted (afs) 
            % field to existing data structures
			for i = 1:numel(sampleData)
				sampleData(i) = subtractAutofluorescence(sampleData(i), autofluor);
			end
			for i = 1:numel(controlData)
				if isempty(controlData(i).(channels{1})), continue, end % % Some FlowData tcData will be empty 
				controlData(i) = subtractAutofluorescence(controlData(i), autofluor);
			end
            
            % Process the compensation - adds matrix compensated (mComp) field
            % to existing data structures
            for i = 1:numel(sampleData)
                sampleData(i) = processCompensation(sampleData(i), coeffs);
            end
            for i = 1:numel(controlData)
				if isempty(controlData(i).(channels{1})), continue, end % % Some FlowData tcData will be empty 
                controlData(i) = processCompensation(controlData(i), coeffs);
            end
            
            % Check compensation results
            [coeffs_comp, ints_comp] = getCoefficients(true);
            
			fitParams = struct( ...
				'intercepts', ints, ...
				'autofluor', autofluor, ...
				'coefficients', coeffs, ...
				'comp_ints', ints_comp, ...
				'comp_coeffs', coeffs_comp);
			
            
            % --- Helper functions --- %
            
            
            function checkInputs_compensateMatrixBatch()
                
                % Check types
                validateattributes(sampleData, {'struct'}, {}, mfilename, 'inputData', 1);
                validateattributes(controlData, {'struct'}, {}, mfilename, 'controlData', 2);
				validateattributes(channels, {'cell', 'char'}, {}, mfilename, 'channels', 3);
				if exist('dataType', 'var')
					validatestring(dataType, fieldnames(controlData(1).(channels{1})), mfilename, 'dataType', 4);
				end
                if exist('gate', 'var')
                    validatestring(gate, fieldnames(controlData(1).gates), mfilename, 'gate', 5);
                end
                if exist('plotsOn', 'var')
                    validateattributes(plotsOn, {'logical'}, {}, mfilename, 'plotsOn', 6);
                else
                    plotsOn = false;
                end
                
                % Check if number of channels is equal to length of single-color struct
                %   (should be 1:1 matching, or fewer if not all are needed)
                if ischar(channels), channels = {channels}; end
                assert((length(channels) + 1) <= length(controlData), ...
                    'Too many channels (%d) for given number of controls (%d)', ...
                    length(channels), length(controlData));
                
                % Check if each given single color channel is in the measured channels
                for ch = 1:numel(channels)
                    if isempty(intersect(channels{ch}, fieldnames(controlData(ch))))
                        error('Given scChannel: %s is not valid', channels{ch});
                    end
                end
            end
            
            
            function autofluor = getAutofluor()
                % Builds an Nx1 array of autofluorescent values in each channel
				%
				%	CURRENTLY NOT USED
                
                % Set up figure to view distribution (if applicable)
                if plotsOn
                    figure()
                    spIdx = 0;
                end
                
                autofluor = zeros(1, numel(channels));
                if exist('gate', 'var')
                    wtGate = controlData(end).gates.(gate);
                else
                    wtGate = true(1, numel(controlData(end).(channels{1}).(dataType)));
                end
                for ch = 1:numel(channels)
                    autofluor(ch) = median(controlData(end).(channels{ch}).(dataType)(wtGate));
                    
                    if plotsOn
                        spIdx = spIdx + 1;
                        ax = subplot(1, numel(channels), spIdx);
                        hold(ax, 'on')
                        histogram(ax, real(controlData(end).(channels{ch}).(dataType)(wtGate)))
                        histogram(ax, real(controlData(end).(channels{ch}).(dataType)(wtGate) - autofluor(ch)))
                    end
                end
            end
            
            
            function dataSubset = subtractAutofluorescence(dataSubset, autofluor)
                % Subtracts autofluorescence from the given data subset
                
                % Iterate over data and subtract autofluorescence
				% Don't gate the data to keep the .afs data full-size.
                for ch = 1:numel(channels)
                    dataSubset.(channels{ch}).afs = ...
							dataSubset.(channels{ch}).(dataType) - autofluor(ch);
                end
            end
            
            
            function [coefficients, intercepts] = getCoefficients(isComp)
                % Finds the coefficients for linear fits between each channel
                % when observing bleed-through.
                %
                % Optional input: isComp (can use to check post-comp data)
                
                % Set up figure to view fitting (if applicable)
                if plotsOn
                    figure()
                    spIdx = 0;
                end
                
                coefficients = zeros(numel(channels));
				intercepts = zeros(numel(channels));
                for chF = 1:numel(channels)
                    
                    % Find regression with each channel
                    for chB = 1:numel(channels)
                        
						if exist('gate', 'var')
							scGate = controlData(chB).gates.(gate);
						else
							scGate = true(size(controlData(chB).(channels{1}).(dataType)));
						end
						
						% Throw out first 10% of cells in case of carry-through
						% --> Cells are ordered by time aquired
% 						timeCutoff = 0.9; % First fraction of samples to ignore
% 						timeGate = true(size(controlData(chB).(channels{1}).(dataType)));
% 						timeGate(1:round(numel(timeGate) * timeCutoff)) = false;
% 						scGate = (scGate & timeGate); 
						
						if (exist('isComp', 'var') && isComp)
							scBleedData = controlData(chB).(channels{chB}).mComp(scGate);
							scFixData = controlData(chB).(channels{chF}).mComp(scGate);
						else
							scBleedData = controlData(chB).(channels{chB}).(dataType)(scGate);
							scFixData = controlData(chB).(channels{chF}).(dataType)(scGate);
						end
						
						% Only necessary for log-scale calculaitons, but linear
						% is the correct way to do it.
% 						valid = (scBleedData > 0 & scFixData > 0);
% 						scBleedData = scBleedData(valid);
% 						scFixData = scFixData(valid);
						
						% The time gating seemed to not get rid of all bad
						% points, causing problems to persist in fitting.
						% Instead, we will try to exclude outliers by egnoring
						% any point that is >10x greater on the fix axis as the
						% bleed axis. 
						badPoints = scFixData > scBleedData .* 10;
						scFixData = scFixData(~badPoints);
						scBleedData = scBleedData(~badPoints);
                        
                        % This is the same as calculating the slope with a least
                        % squares metric. Rows are channels being bled into and
                        % columns are bleeding channels. 
%                         coefficients(chF, chB) = scBleedData(:) \ scFixData(:);

						% Best to also gather the intercept so we can estimate
						% the autofluorescence better
                        [~, coefficients(chF, chB), intercepts(chF, chB)] = regression(scBleedData', scFixData');
                        
                        if plotsOn
							maxPoints = 6000;
							numPoints = min(numel(scFixData), maxPoints);
							points = randperm(numel(scFixData), numPoints);
							switch dataType
								case {'mef', 'mefl'}
									
									xrange = logspace(-1, 9, 100);
									
									spIdx = spIdx + 1;
									ax = subplot(numel(channels), numel(channels), spIdx);
									hold(ax, 'on')
									
									plot(ax, Transforms.lin2logicleMEF(scBleedData(points)), ...
										 Transforms.lin2logicleMEF(scFixData(points)), ...
										 '.', 'MarkerSize', 4)
									plot(ax, Transforms.lin2logicleMEF(xrange), ...
										 Transforms.lin2logicleMEF(xrange * coefficients(chF, chB) + intercepts(chF, chB)), ...
										 '-', 'linewidth', 4)
									Plotting.biexpAxesMEF(ax);
																	
								otherwise
									
									xrange = logspace(-1, 6, 100);
									
									spIdx = spIdx + 1;
									ax = subplot(numel(channels), numel(channels), spIdx);
									hold(ax, 'on')

									plot(ax, Transforms.lin2logicle(scBleedData(points)), ...
										 Transforms.lin2logicle(scFixData(points)), ...
										 '.', 'MarkerSize', 4)
									plot(ax, Transforms.lin2logicle(xrange), ...
										 Transforms.lin2logicle(xrange * coefficients(chF, chB) + intercepts(chF, chB)), ...
										 '-', 'LineWidth', 4)
									Plotting.biexpAxes(ax);
							end
							
							% Axis labeling
							title(sprintf('Slope: %.2f | Intercept: %.2f', ...
								coefficients(chF, chB), intercepts(chF, chB)), 'fontsize', 14)
							if (chF == numel(channels))
								xlabel(strrep(channels{chB}, '_', '-'))
							end
							if (chB == 1)
								ylabel(strrep(channels{chF}, '_', '-'))
							end
                        end
                    end
                end
            end

            
            function dataSubset = processCompensation(dataSubset, coefficients)
                % Processes the compensation for one data subset
                
                % Extract raw data into NxM matrix (N = num channels, M = num cells)
				% Don't gate the data to keep the .mComp data full-size.
				rawData = zeros(numel(channels), numel(dataSubset.(channels{1}).afs));
                for ch = 1:numel(channels)
					rawData(ch, :) = dataSubset.(channels{ch}).afs;
                end
                
                % Apply matrix compensation
                fixedData = Compensation.compensateMatrix(rawData, coefficients);
                
                % Implant compensated data
                for ch = 1:numel(channels)
                    dataSubset.(channels{ch}).mComp = fixedData(ch, :)';
                end
                
%                 % Test plotting
%                 if exist('plotsOn','var') && plotsOn
%                     figure()
%                     biexplot(scBleedChannel, scFixChannel, {'.', 'markersize', 5})
%                     hold on
%                     biexplot(inBleedChannel, inFixChannel, {'.', 'markersize', 5})
%                     hold on
%                     biexplot(inBleedChannel, compdInFixChannel, {'.', 'markersize', 5})
%                     pause()
%                 end

            end
        end
            
        
        function compdData = compensateMatrix(inputData, coefficients)
            % compensateMatrix compensates data using linear coefficients for
            % bleed-trhough between channels.
            % 
            %   compData = compensateMatrix(matrixFile, inputData)
            %   
            %       Inputs: 
            %           inputData       An NxM matrix of data points corresponding with
            %                           M cells in N channels to be compensated. 
            %
            %           coefficients    An NxN matrix of coefficients corresponding with
            %                           the slope of lines of best fit between N channels 
            %                           bleeding into each other. 
            %
            %       Outputs:
            %           compData        A copy of the inputData struct with added matComp field
            %            (struct)       for each compensated channel.
            %   
            %   This method solves the linear set of equations:
            %       X_observed = A * X_real
            %
            %   'inputData' corresponds with 'X_observed' and 'coefficients'
            %   corresponds with the matrix 'A'. We invert 'A' to solve for
            %   X_real, which is returned as 'compdData'.
            %
            % Ross Jones
            % Weiss Lab
            % 2017-06-06
            %
            % Update Log:
            %
            % ...
            
            % Check inputs are valid
            checkInputs_compensateMatrix(inputData, coefficients);
            
            % Compensate data
            compdData = coefficients \ inputData;
            
            
            % --- Helper functions --- %
            
            
            function checkInputs_compensateMatrix(inputData, coefficients)
                
                % Check types
                validateattributes(inputData, {'numeric'}, {}, mfilename, 'inputData', 1);
                validateattributes(coefficients, {'numeric'}, {}, mfilename, 'coefficients', 2);
                
                % Check sizes
                assert(size(inputData, 1) == size(coefficients, 2), ...
                    '# Channels in ''inputData'' (%d) and ''coefficients'' (%d) are not the same!', ...
                    size(inputData, 1), size(coefficients, 2))
                assert(size(coefficients, 1) == size(coefficients, 2), ...
                    '''coefficients'' is not square! (H = %d, W = %d)', ...
                    size(coefficients, 1), size(coefficients, 2))
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