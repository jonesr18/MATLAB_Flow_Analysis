classdef FlowAnalysis < handle
    % Compilation of analytical procedures from the Weiss Lab flow cytometry MATLAB repository
	%
	%	Methods are implemented as static functions, so they can be called directly
	%	or after creating a handle to this class.
	%
	%		Ex:
	%		data = FlowAnalysis.openFiles(TF_marker, args);
	%
	%		Ex2:
	%		flow = FlowAnalysis();
	%		data = flow.bin(data, const, out);
	%
    %   Functions:
    %
    %       data = openFiles(varargin)
    %       [xBleed, yComp] = compensateSC(JCFile, singleColorFile, correctFile, bC, cC, plotsOn)
    %       compData = compensateBatchSC(jcData, scData, channels, inputData, gate, plotsOn)
    %       compdData = compensateMatrix(matrixFile, inputData)
    %       YI = lsq_lut_piecewise(x, y, XI)
    %       data = bin(data, constitutive, outputs, gate)
    %       data = cluster(data, channels, method, dataType, numPop)
    %       [data, stats] = calcStats(inputs [data, channels, dataType, threshChan, thresh, threshPct, statTypes] )
    %       data = mixtureModel(data, channels, dataType, numPop, par)
    %       [x1, x2, p1, p2, yfit1, yfit2, rsq_both] = LSB_fit(constitutive, reporter, doLog, skip)
    %       data = ronBin(data, constitutive, outputs)
    %
    % Written/Compiled by Ross Jones
    % Weiss Lab, MIT
    % Last updated 2016-05-27
    
    methods (Static)
        
        function data = openFiles(varargin)
            % Opens the given files (opens a GUI to select files if no args)
            %
            % The following is no longer true, but kept in case it relevant portions are uncommented
            %   The first argument is the name of the channel with the transfection
            %   marker. Note that channel names are changed such that ' ' and '-'
            %   become '_', and '#' becomes 'N', so take that into consideration
            %   when passing the TF_marker channel name. This is used to sort the
            %   data based on the 'raw' .fcs data.
            %
            % Inputs:
            %   Subsequent arguments must be valid filenames given as strings.
            %   .fcs at the end of the filenames can be left off. 
            %   If files are on path, name can be given, otherwise full path + name
            %   is required. An error will occur if this is not followed.
            %
            % 
            % Output:
            %   A standard struct of data with the following characteristics:
            %       Fields:
            %           header          The header copied from the .fcs files
            %           nObs            The number of observations (events)
            %           chanNames       A list of channels, given as strings. 
            %                               Channel names are optimized for use as fields, see below
            %           <channelName>   Each channel name in chanNames is also a field, pointing 
            %                           to a substructure. The substructure has the fields 'raw', 
            %                           'scaled', and 'comp'. Comp is only useful if you did 
            %                           compensation on the machine. Scaled is also by the machine.
            %                           If comp and scaled are the same as raw, they are not
            %                           included in the data strcuture. 
            %
            %   Relevant header fields:
            %       filename
            %       TotalEvents
            %       CompLabels
            %       CompMat
            %       par.name         (parameters -- channels)
            %
            %
            % Written by:
            %   Ross Jones
            %   Weiss Lab, MIT
            %   Updated 2015-11-18
            
            % Check inputs (should be the transfection marker channel name, then filenames)
%             validateattributes(TF_marker, {'char'}, {}, mfilename, 'TF_marker', 1);
            filenames = cell(1, nargin);
            for i = 1: numel(varargin)
                name = varargin{i};
                if (~strcmpi(name(end-3:end), '.fcs'))
                    name = strcat(name, '.fcs');
                end
                if (~exist(name, 'file'))
                    error('Filename %s does not exist in the path', name);
                end
                filenames{i} = name;
            end
            
            % If no filenames were given, open GUI to find files
            if (nargin == 0)
                [filenames, filepath] = uigetfile('*.fcs', 'Select FCS file', 'MultiSelect', 'on');
            end
            
            % Data will be stored in a handy struct
            data = struct();

            % Check if multiple files were selected, if not, put the lone filename
            % into a cell array so that it can be processed the same.
            if ~iscell(filenames)
                filenames = {filenames};
            end
            numFiles = numel(filenames);
            
            % Run through and extract data from files, putting it into the
            % standardized struct array.
            for i = 1:numFiles
                % Depending on how the files were added, append the filepaths
                if (nargin == 0)
                    % In this case, they were selected from the GUI, so we need the path
                    fullfile = strcat(filepath, filenames{i});
                else
                    % In this case, they were given, so we assume they are on the path
                    % or the full path was given with the file.
                    fullfile = filenames{i};
                end
                [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(fullfile);

                % Get channel names and store relevant data in data's fields
                channelNames = {fcshdr.par.name};
                data(i).header = fcshdr;
                data(i).nObs = size(fcsdat, 1);
                
                % Ensure channel names are legal field declarations
                for channel = 1:numel(channelNames)
                    chanName = channelNames{channel};
                    chanName(chanName == '-') = '_';
                    chanName(chanName == ' ') = '_';
                    chanName(chanName == '#') = 'N';
                    channelNames{channel} = chanName;
                end
                
                % I've erased some fucntionality that I used to have where everything 
                % was sorted based on a transfection marker, which is unnecessary.
                
%                 % Check that TF_marker exists
%                 TF_marker_idx = find(strcmpi(TF_marker, channelNames), 1);
%                 if (isempty(TF_marker_idx))
%                     error('TF_marker not found in file.\nTF_marker: %s\n', TF_marker);
%                 end
%                 
%                 % Sort based on the values in the TF_marker
%                 [~, sortIDX] = sort(fcsdat(:, TF_marker_idx));
                
                % Extract data :: The data is in columns, so we need the column numbers
                for channel = 1:numel(channelNames)
                    
                    chanName = channelNames{channel};
                    ch = struct();
                    
                    % Extract raw data
%                     ch.raw = fcsdat(sortIDX, channel);
                    ch.raw = fcsdat(:, channel);
                    
                    % Only taken if the scaled is different
                    if (~isempty(fcsdatscaled) && all(fcsdatscaled(:, channel) ~= fcsdat(:, channel)))
    %                     ch.scaled = fcsdatscaled(sortIDX, channel);
                        ch.scaled = fcsdatscaled(:, channel);
                    end
                    
                    % only taken if the comp is different
                    if (~isempty(fcsdatcomp) && all(fcsdatcomp(:, channel) ~= fcsdat(:, channel)))
%                         ch.comp = fcsdatcomp(sortIDX, channel);
                        ch.comp = fcsdatcomp(:, channel);
                    end
                    
                    % Store channel data
                    data(i).(chanName) = ch;
                end

                % Has adjusted names
                data(i).chanNames = channelNames;
            end

            fprintf(1, 'Finished creating data struct\n');
        end
		
		
		function [xBleed, yComp] = compensateSC(JCFile, singleColorFile, correctFile, bC, cC, plotsOn)
			% [xBleed,yComp]=compensate(JCFile, singleColorFile, correctFile, bC, cC,plotsOn)
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
			YI = FlowAnalysis.lsq_lut_piecewise(bValsSC, cValsSC, XI);
			
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
        
        
        function compData = compensateBatchSC(jcData, scData, scChannels, inputData, fixChannels, gate, plotsOn)
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
            parseInputs();
                        
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
            
            
            function parseInputs()
                
                % Check types
                validateattributes(jcData, {'struct'}, {}, mfilename, 'jcData', 1);
                validateattributes(scData, {'struct'}, {}, mfilename, 'scData', 2);
                validateattributes(scChannels, {'cell', 'char'}, {}, mfilename, 'channels', 3);
                validateattributes(inputData, {'struct'}, {}, mfilename, 'inputData', 4);
                validateattributes(fixChannels, {'cell', 'char'}, {}, mfilename, 'channels', 5);
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
                YI = FlowAnalysis.lsq_lut_piecewise(scBleedChannel, scFixChannel, XI);
                
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
        
        
        function [data, jcData, scData, autofluor, coeffs, coeffs_comp] = ...
                        compensateMatrixBatch(data, channels, jcData, scData, gate, plotsOn)
            % compensateMatrixBatch compensates an entire standard data struct by generating 
            % coefficients with linear fitting on single-color control data and processing the
            % data to be sent piecewise to FlowAnalysis.compensateMatrix(). 
            % 
            % The method first subtracts autofluorescence based on the untransfected cell
            % sample provided (see below)
            %
            %   compData = compensateMatrixBatch(justCellsFile, singleColorFile, colors, inputData)
            %   
            %       Inputs: 
            %           data            A stardard struct with channel data to be compensated. 
            %            (struct)       'fixChannels' should all exist in this struct.
            %           
            %           channels        A cell array of strings indicating the channels to be fixed
            %            (cell)         by the compensation procedure.
            %
            %                           Entries should be strings with the names given by Fortessa
            %                           but with spaces and dashes replaced by underscores
            %                               eg: {'PE_Texas_Red_A', 'FITC_A'}
            %                               Note: this is case-sensitive
            %
            %           jcData          A standard struct with untransfected cells data
            %            (struct)       
            %
            %           scData          A standard struct with single color samples data
            %            (struct)       
            %                       *** The order MUST BE THE SAME as the corresponding 
            %                           colors in 'channels'. This will not be auto-checked!
            %
            %           gate (opt)      A string indicating what gated population to compensate
            %            (string)       with, for example, 'P3'. Must match entries in 'gates'
            %                           field of all data structs.
            %
            %           plotsOn (opt)   A boolean indicating whether to plot or not
            %            (logical)
            %
            %       Outputs:
            %           data            A copy of the data struct with added 'afs' and 'mComp' 
            %            (struct)       fields for each compensated channel.
            %                               afs = autofluorescence subtracted
            %                               mComp = matrix compensated
            %
            %           jcData          A copy of the jcData struct with added 'afs' and 'mComp' 
            %            (struct)       fields for each compensated channel.
            %
            %           scData          A copy of the scData struct with added 'afs' and 'mComp' 
            %            (struct)       fields for each compensated channel.
            %
            %
            % Ross Jones, Weiss Lab
            % jonesr18@mit.edu
            % Created 2017-06-06
            % 
            % Update Log: 
            %
            % ...
            
            % Check the inputs - makes no changes to data structures
            checkInputs();
            
            % Subtract autofluorescence - adds autofluorescence subtracted (afs) 
            % field to existing data structures
            
            % Find medians for each channel in jcData
            autofluor = getAutofluor();
            for i = 1:numel(data)
                data(i) = subtractAutofluorescence(data(i), autofluor);
            end
            for i = 1:numel(scData)
                scData(i) = subtractAutofluorescence(scData(i), autofluor);
            end
            jcData = subtractAutofluorescence(jcData, autofluor);
            
            % Find coefficients using least-squares linear fits - uses afs data
            coeffs = getCoefficients();
            
            % Process the compensation - adds matrix compensated (mComp) field
            % to existing data structures
            for i = 1:numel(data)
                data(i) = processCompensation(data(i), coeffs);
            end
            for i = 1:numel(scData)
                scData(i) = processCompensation(scData(i), coeffs);
            end
            jcData = processCompensation(jcData, coeffs);
            
            % Check compensation results
            coeffs_comp = getCoefficients(true);
            
            
            % --- Helper functions --- %
            
            
            function checkInputs()
                
                % Check types
                validateattributes(data, {'struct'}, {}, mfilename, 'inputData', 1);
                validateattributes(channels, {'cell', 'char'}, {}, mfilename, 'channels', 2);
                validateattributes(jcData, {'struct'}, {}, mfilename, 'jcData', 3);
                validateattributes(scData, {'struct'}, {}, mfilename, 'scData', 4);
                if exist('gate', 'var')
                    validatestring(gate, fieldnames(jcData(1).gates), mfilename, 'gate', 5);
                end
                if exist('plotsOn', 'var')
                    validateattributes(plotsOn, {'logical'}, {}, mfilename, 'plotsOn', 6);
                else
                    plotsOn = false;
                end
                
                % Check if number of channels is equal to length of single-color struct
                %   (should be 1:1 matching, or fewer if not all are needed)
                if ischar(channels), channels = {channels}; end
                assert(length(channels) == length(scData), ...
                    'Number of channels (%d) does not equal the number of single-color controls (%d)', ...
                    length(channels), length(scData));
                
                % Check if each given single color channel is in the measured channels
                for ch = 1:numel(channels)
                    channel = channels{ch};
                    if isempty(intersect(channel, fieldnames(scData(ch))))
                        error('Given scChannel: %s is not valid', channel{:});
                    end
                end
            end
            
            
            function autofluor = getAutofluor()
                % Builds an Nx1 array of autofluorescent values in each channel
                
                % Set up figure to view distribution (if applicable)
                if (exist('plotsOn', 'var') && plotsOn)
                    figure()
                    spIdx = 0;
                end
                
                autofluor = zeros(1, numel(channels));
                if exist('gate', 'var')
                    jcGate = jcData.gates.(gate);
                else
                    jcGate = true(1, numel(jcData.(channels{1}).raw));
                end
                for ch = 1:numel(channels)
                    autofluor(ch) = median(jcData.(channels{ch}).raw(jcGate));
                    
                    if (exist('plotsOn', 'var') && plotsOn)
                        spIdx = spIdx + 1;
                        ax = subplot(1, numel(channels), spIdx);
                        hold(ax, 'on')
                        histogram(ax, jcData.(channels{ch}).raw(jcGate))
                        histogram(ax, jcData.(channels{ch}).raw(jcGate) - autofluor(ch))
                    end
                end
            end
            
            
            function dataSubset = subtractAutofluorescence(dataSubset, autofluor)
                % Subtracts autofluorescence from the given data subset
                
                % Iterate over data and subtract autofluorescence
                if exist('gate', 'var')
                    dGate = dataSubset.gates.(gate);
                else
                    dGate = true(1, numel(dataSubset.(channels{1}).raw));
                end
                for ch = 1:numel(channels)
                    dataSubset.(channels{ch}).afs = dataSubset.(channels{ch}).raw(dGate) - autofluor(ch);
                end
            end
            
            
            function coefficients = getCoefficients(isComp)
                % Finds the coefficients for linear fits between each channel
                % when observing bleed-through.
                %
                % Optional input: isComp (can use to check post-comp data)
                
                % Set up figure to view fitting (if applicable)
                if (exist('plotsOn', 'var') && plotsOn)
                    figure()
                    spIdx = 0;
                end
                
                coefficients = zeros(numel(channels));
                for chF = 1:numel(channels)
                    
                    % Find regression with each channel
                    for chB = 1:numel(channels)
                        
                        if (exist('isComp', 'var') && isComp)
                            scBleedData = scData(chB).(channels{chB}).mComp;
                            scFixData = scData(chB).(channels{chF}).mComp;
                        else
                            scBleedData = scData(chB).(channels{chB}).afs;
                            scFixData = scData(chB).(channels{chF}).afs;
                        end
                        
                        % This is the same as calculating the slope with a least
                        % squares metric. Rows are channels being bled into and
                        % columns are bleeding channels. 
%                         coefficients(chF, chB) = scBleedData(:) \ scFixData(:);
                        [~, coefficients(chF, chB)] = regression(scBleedData', scFixData');
                        
                        if (exist('plotsOn', 'var') && plotsOn)
                            
                            subSample = linspace(min(scBleedData), max(scBleedData), 100);
                            
                            spIdx = spIdx + 1;
                            ax = subplot(numel(channels), numel(channels), spIdx);
                            hold(ax, 'on')
                            
                            plot(ax, Transforms.lin2logicle(scBleedData(1:10:end)), ...
                                Transforms.lin2logicle(scFixData(1:10:end)), ...
                                '.', 'MarkerSize', 2)
                            
                            plot(ax, Transforms.lin2logicle(subSample), ...
                                Transforms.lin2logicle(subSample * coefficients(chF, chB)), ...
                                '-', 'linewidth', 4)
                            Plotting.biexpAxes(ax);
                            if (chB ~= chF)
                                ax.YLim = Transforms.lin2logicle([-100, 200]);
                            end
                        end
                    end
                end
            end

            
            function dataSubset = processCompensation(dataSubset, coefficients)
                % Processes the compensation for one data subset
                
                % Extract raw data into NxM matrix (N = num channels, M = num cells)
                rawData = zeros(numel(channels), numel(dataSubset.(channels{1}).afs));
                for ch = 1:numel(channels)
                    rawData(ch, :) = dataSubset.(channels{ch}).afs;
                end
                
                % Apply matrix compensation
                fixedData = FlowAnalysis.compensateMatrix(rawData, coefficients);
                
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
            checkInputs(inputData, coefficients);
            
            % Compensate data
            compdData = coefficients \ inputData;
            
            
            % --- Helper functions --- %
            
            
            function checkInputs(inputData, coefficients)
                
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
        
        
        function data = bin(data, inputs, outputs, dataType, numBins, binMin, binMax)
            % data = bin(data, constitutive, outputs)
            %
            % Creates bins for the given data, adding the 'bins' and 'gmeans' fields
            % to the data struct. 
            %
            %   Inputs
            %       'data'          A data struct like that generated by the openFiles() method.
            %       'inputs'        A cell array of string input channel names. Must match 
            %                       a field in the data struct with the subfield 'raw'.
            %                       For single channel inputs, a string can be given.
            %       'outputs'       Similar to inputs, but for output channels. Can be a cell 
            %                       array of channel names or a single channel name string.
            %       'dataType'      The data type to use, eg 'raw', or 'scComp'
            %       'numBins'       (Optional) The number of bins to use - same for each dimension. 
            %                           < Default: 20 >
            %       'binMin'        (Optional) A double value indicating the
            %                       lower edge of the smallest bin 
            %                           < Default: lin2logicle(-171.088) >
            %       'binMax'        (Optional) A double value indicating the
            %                       upper edge of the largest bin
            %                           < Default: lin2logicle(2^18) >
            %
            %   Outputs
            %       'data.bins'             A boolean matrix where each row corresponds
            %                               with a bin and each column is a cell. An
            %                               element is true if the cell is in the bin.
            %       'data.(...).gmeans'     Bin geometric means (positive values only)
            %       'data.(...).medians'    Bin median
            %       'data.(...).stds'       Bin std deviation
            %       'data.(...).geostds'    Bin geometric standard deviation (pos vals only)
            %
            % Updated 2017-05-15
            % - Changed data.bins to be a double array of bin IDs instead of a cell
            %   array of logicle vectors
            
            % Check inputs
            checkInputs()
            
            % Collect all channels
            channels = [inputs, outputs];
            
            % Find binning dimensions
            nDims = numel(inputs);
            binDims = numBins * ones(1, nDims);
            if nDims == 1
                binDims = [1, binDims];
            end
            
            % Create bins with input channels
            for i = 1:numel(data)
                data(i).bins = generateBins(data(i), nDims);
            end
            
            fprintf(1, 'Finished creating bins\n');
            
            % Calculate bin statistics
            for i = 1:numel(data)
                
                % Pre-acllocate substructure arrays for gmeans calculation
                for chan = channels
                    data(i).(chan{:}).gmeans = nan(binDims);
                    data(i).(chan{:}).medians = nan(binDims);
                    data(i).(chan{:}).stds = nan(binDims);
                    data(i).(chan{:}).geostds = nan(binDims);
                end

                % Iterate over bins, finding geomeans of bins 
                for j = 1:numel(data(i).bins)
                    bin = data(i).bins{j};
                    
                    % Check to make sure bin isn't empty
                    if isempty(bin)
                        continue
                    end
                    
                    % Since there are a variable number of outputs, we have to 
                    % handle them all and ensure no negative values are being
                    % included, as geometric means don't work for negatives!
                    nonNeg = true(numel(bin), 1);
                    chanData = cell(1, numel(channels));
                    for ch = 1:numel(channels)
                        chanData{ch} = data(i).(channels{ch}).(dataType)(bin);
                        nonNeg = (nonNeg & (chanData{ch} > 0));
                    end

                    % Calculate channel metrics
                    for ch = 1:numel(channels)
                        data(i).(channels{ch}).medians(j) = median(chanData{ch}, 'omitnan');
                        data(i).(channels{ch}).stds(j) = std(chanData{ch}, 'omitnan');

                        if (sum(nonNeg) > 0) % Handle if there are non non-neg vals
                            data(i).(channels{ch}).gmeans(j) = geomean(chanData{ch}(nonNeg));
                            data(i).(channels{ch}).geostds(j) = geostd(chanData{ch}(nonNeg));
                        else
                            data(i).(channels{ch}).gmeans(j) = nan;
                            data(i).(channels{ch}).geostds(j) = nan;
                        end
                    end
                end
            end
            
            
            % --- Helper Functions --- %
            
            function checkInputs()
                
                validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                validateattributes(inputs, {'char', 'cell'}, {}, mfilename, 'inputs', 2);
                validateattributes(outputs, {'char', 'cell'}, {}, mfilename, 'outputs', 3);
                if (~iscell(inputs))
                    inputs = {inputs};      % This just makes things less messy throughout
                end
                if (~iscell(outputs))
                    outputs = {outputs};    % This just makes things less messy throughout
                end
                validatestring(dataType, fieldnames(data(1).(inputs{1})), mfilename, 'dataType', 4);

                % Check optional parameters
                if ~exist('numBins', 'var')
                    numBins = 20;   % Default 20 bins
                end
                if ~exist('binMin', 'var')
                    binMin = 0;     % logicle2lin(0) = -171.088
                end
                if ~exist('binMax', 'var')
                    binMax = 4.5;   % logicle2lin(4.5) = 2^18
                end
            end
            
            
            function bins = generateBins(dataToBin, nDims)
                % Generates bins for the given data sample in N dimensions where
                % N is given by nDims
                
                % Extract cells
                cells = zeros(numel(dataToBin.(inputs{1}).(dataType)), nDims);
                for d = 1:nDims
                    cells(:, d) = Transforms.lin2logicle(dataToBin.(inputs{d}).(dataType));
                end
                
                % Do binning
                bins = FlowAnalysis.simpleBin(cells, numBins, binMin, binMax);
            end
        end
        
        
        function bins = simpleBin(cells, numBins, binMin, binMax)
            % Takes a double array of cells and bins them in biexp space, 
            % returning a double array of bin IDs for each cell
            %
            %   Inputs
            %       'cells'         A N x M double matrix of N cell fluorescence  
            %                       values in N pre-selected channels. If using 
            %                       default parameters, this matrix should be 
            %                       logicle-transformed.
            %       'numBins'       (Optional) A double value indicating the
            %                       number of bins to use 
            %                           < Default: 20 >
            %       'binMin'        (Optional) A double value indicating the
            %                       lower edge of the smallest bin 
            %                           < Default: lin2logicle(-171.088) >
            %       'binMax'        (Optional) A double value indicating the
            %                       upper edge of the largest bin
            %                           < Default: lin2logicle(2^18) >
            %
            %   Output
            %       'bins'          A N-dimensional cell matrix where each cell 
            %                       holds a double array of numerical indexes 
            %                       for each cell in a given bin. The size of
            %                       each dimension is given by 'numBins'.
            
            % Check optional parameters
            if ~exist('numBins', 'var')
                numBins = 20;   % Default 20 bins
            end
            if ~exist('binMin', 'var')
                binMin = 0;     % logicle2lin(0) = -171.088
            end
            if ~exist('binMax', 'var')
                binMax = 4.5;   % logicle2lin(4.5) = 2^18
            end
            
            % Define bins in biexponential space at a wide range
            edges = linspace(binMin, binMax, numBins + 1);
            
            % Setup bins as a cell array with dimensions equal to number of cell channels. 
            nDims = size(cells, 2);
            binDims = numBins * ones(1, nDims);
            if nDims == 1
                binDims = [binDims, 1];
            end
            bins = cell(binDims);
            binsNum = reshape(1:numel(bins), binDims);
            
            % Find cells in each bin
            for i = 1:numel(bins)
                
                % Find bin coordinates in ND-space
                [I, J] = find(binsNum == i);
                while (J > numBins)
                    [I2, J] = find(binsNum == J);
                    I = [I, I2]; %#ok<AGROW>
                end
                binCoords = [I, J];
                
                % Identify all cells within bin by checking bin edges in each dimension
                cellsInBin = true(size(cells(:, 1)));
                for j = 1:nDims
                    cellsInBin = (cellsInBin & ...
                                 (cells(:, j) <= edges(binCoords(j) + 1)) & ...
                                 (cells(:, j) > edges(binCoords(j))));
                end
                bins{i} = find(cellsInBin);
            end
            
        end
        
        
        function data = cluster(data, channels, method, dataType, numPop)
            % CLUSTER Clusters the given data with the given clustering method
            %
            %   Parameters:
            %       
            %       data        Struct: Standard data structure
            %       channels    Cell<String>: A list of names of channels to fit on
            %                    - will also accept a single string for evaluating one channel
            %       method      String: The clustering method to use
            %                    - Valid entities: 'kmeans', 'svm', 'otsu', 'gmm'
            %                    - Can be Cell<String> for doing multiple clustering methods!
            %                    - kmeans, svm, and gmm will use all channels simultaneously to 
            %                      compute clusters, whereas otsu will use each channel individually
            %                    - To force channels to be assessed independently, run this method
            %                      sequentially with different channels given individually.
            %       dataType    (optional) String: The type of data to fit 
            %                    - default = 'raw'
            %       numPop      (optional) Int: The number of populations to fit 
            %                    - For k-means and gmm only! 
            %                    - default = 2
            %
            %   Ouputs:
            %       data        Updated data struct. Each fitted channel now has a field with the
            %                   same name as the one given for the 'method' parameter with further
            %                   sub-fields with information as such:
            %
            %                   kmeans:
            %                       idx         (indexes of classes)
            %                       centroids   (centroids of classes)
            %                   svm:
            %                       <not implemented>
            %                       
            %                   otsu:
            %                       thresh      (threshold value)
            %                       pos         (positive cells)
            %                       neg         (negative cells)
            %                       em          (effectiveness matrix)
            %                   gmm:
            %                       mus         (means of populations)
            %                       sigmas      (covariance of populations)
            %                       idx         (indexes of classes)
            %                       post        (posterior probability of being in each class)
            %                       pdf         (probability density function)
            %                                   - eval at linspace(-2, 5, 1e2) for each channel
            %    
            % Written by Ross Jones
            % jonesr18@mit.edu
            % Weiss Lab, MIT
            % Last updated 2016-04-05
            
            % Check existence of optional inputs, assign defaults
            if (~exist('dataType', 'var'))
                dataType = 'raw';
            end
            if (~exist('numPop', 'var'))
                numPop = 2;
            end
            
            % Check inputs
            if (ischar(channels))
                channels = {channels};
            end
            validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
            validateattributes(channels, {'cell'}, {'vector'}, mfilename, 'channels', 2);
            s = setdiff(channels, fieldnames(data(1)));
            if (~isempty(s))
                error('Invalid channel given: %s\n', s{:});
            end
            s = setdiff(method, {'kmeans', 'svm', 'otsu', 'gmm'});
            if (~isempty(s))
                error('Invalid method given: %s\n', s{:});
            end
            validatestring(dataType, fieldnames(data(1).(channels{1})), mfilename, 'dataType', 4);
            validateattributes(numPop, {'numeric'}, {'positive'}, mfilename, 'numPop', 5);
            numPop = round(numPop);
            
            % Find data sizes
            D = numel(data);
            C = numel(channels);
            
            % For tracking progress
            totalProgress = D;
            currProgress = 0;
            percentile = 10;
            
            
            % Set fitting parameters
            options = statset('MaxIter', 500);
            kmeansParams = {
                'Distance', 'sqeuclidean', ...  % Use squared-euclidean distance metric
                'Start', 'plus', ...            % Use k-means++ to seed initial guesses
                'Replicates', 10};              % 10 replicates to avoid local minima
            gmmParams = {
                'Start', 'plus', ...            % Use k-means++ to seed initial guesses
                'Replicates', 10, ...           % 10 replicates to avoid local minima
                'CovarianceType', 'full', ...   % Covariance of each component is not independent
                'SharedCovariance', false};     % Populations do not share covariance
            
            
            % Automatically create a C-dimensional coordinate grid to deliver pdf points for gmm
            pdfRange = linspace(-2, 5, 1e2);
            if (C == 1)
                coords = pdfRange;
            else
                eval(strcat( '[x1', sprintf(', x%d', 2:C), '] = ndgrid(pdfRange);' ));
                eval(strcat( 'coords = [reshape(x1, [], 1, 1)', ...
                             sprintf(', reshape(x%d, [], 1, 1)', 2:C), '];' ));
            end
            
            
            tf = Transforms();
            for i = 1:D
                
                % Get the total number of cells
                N = numel(data(i).(channels{1}).(dataType));    
                
                % Do different clustering methods
                
                % --- Otsu's method for thresholding ---% 

                if (any(strcmpi(method, 'otsu')))
                    % Otsu is done differently than the others, since it can only look at one
                    % channel at a time.

                    for j = 1:C
                        dataArray = tf.lin2logicle(data(i).(channels{j}).(dataType));
                        edges = linspace(-5, 5, 26);
                        [thresh, effectiveness] = otsuthresh(histcounts(dataArray, edges));
                        adjThresh = thresh * max(dataArray);
                        data(i).(channels{j}).otsu.thresh = adjThresh;
                        data(i).(channels{j}).otsu.pos = dataArray(dataArray >= adjThresh);
                        data(i).(channels{j}).otsu.neg = dataArray(dataArray < adjThresh);
                        data(i).(channels{j}).otsu.em = effectiveness;
                    end
                end
                
                % --- Setup dataArray --- %
                
                if (any(strcmpi(method, {'kmeans', 'svm', 'gmm'})))
                    % The other arrays can be handled normally
                    
                    % Pre-allocate an array with N rows and C columns, where N is the number
                    % of cells in the population and C is the number of channels to look at
                    % - Log transform data first as well
                    dataArray = zeros(N, C);
                    for j = 1:C
%                         fprintf('%d %d\n', i, j);
                        dataArray(:, j) = tf.lin2logicle(data(i).(channels{j}).(dataType));
                    end
                end

                % Use strcmpi instead of switch/case since method can be Cell<String>

                % --- K-means ---% 

                if (any(strcmpi(method, 'kmeans')))

                    % Do k-means clustering
                    [idx, cen] = kmeans(dataArray, numPop, 'Options', options, kmeansParams{:});

                    % Store data
                    for j = 1:C
                        data(i).(channels{j}).kmeans.idx = idx;
                        data(i).(channels{j}).kmeans.centroids = cen(:, j);
                    end
                end

                % --- Support Vector Machine --- %

                if (any(strcmpi(method, 'svm')))
                    warning('SVM not implemented!');
                end

                % --- Gaussian mixture model --- %

                if (any(strcmpi(method, 'gmm')))

                    % Fit Gaussian mixture model
                    gmm = fitgmdist(dataArray, numPop, 'Options', options, gmmParams{:});

                    % Collect pdfs
                    if (C == 1)
                        data(i).(channels{j}).gmm.pdf = gmm.pdf(coords);
                    else
                        % This comes out in a weird single-column format
                        pdfRaw = gmm.pdf(coords);

                        % Fix it by reshaping into a meshgrid the size of the pdfRange^C
                        pdfReshaped = reshape(pdfRaw, repmat(numel(pdfRange), [1, C]));

                        % Permute the X- and Y-dimensions (1 and 2) because for some reason
                        % meshgrids have X changing along the 2nd dimension, which is
                        % undesirable for summing across all other dimensions.
                        pdfPermuted = permute(pdfReshaped, [2, 1, 3:C]);

                        % Sum over all non-j dimensions to get the PDF of a single channel
                        for j = 1:C
                            pdfSummed = pdfPermuted;
                            for k = 1:C
                                if (k ~= j)
                                    pdfSummed = sum(pdfSummed, k);
                                end
                                data(i).(channels{j}).gmm.pdf = pdfSummed;
                            end
                        end
                    end

                    % Collect mus / covars / clusters
                    idx = gmm.cluster(dataArray);
                    post = gmm.posterior(dataArray);
                    for j = 1:C
                        data(i).(channels{j}).gmm.mus = gmm.mu(:, j);
                        data(i).(channels{j}).gmm.sigmas = gmm.Sigma(:, j, :);
                        data(i).(channels{j}).gmm.idx = idx;
                        data(i).(channels{j}).gmm.post = post;
                    end
                end
                % Update progress
                currProgress = currProgress + 1;
                pctProgress = currProgress / totalProgress * 100;
                if (pctProgress > percentile)
                    fprintf(1, '%d%% Complete\n', percentile);
                    percentile = percentile + 10;
                end
            end
        end
        
        
        function [data, stats] = calcStats(inputs)
            % Calculates population statistics for the given data
            %
            %   Inputs (struct array w/ fields below):
            %
            %       data            A standard struct with the data to be analyzed
            %                        - Must have 'channel' and 'dataType' given by user
            %       channels        The channel(s) to calculate stats for
            %       dataType        The data type to use, eg 'raw', or 'scComp' 
            %       threshChan      (optional) The channel(s) to threshold/gate on. 
            %                        - Can be a string or a cell array of strings
            %       thresh          (optional) The value(s) to use for thresholding
            %                        - If more than one threshChan, can give an array of thresh vals
            %       threshPct       (optional) Threshold based on average number of cells above the
            %                       given thresholding value for each sample (for normalization). 
            %       statTypes       (optional) Only returns the indicated fields in 'stats' output
            %                        - Can be a single string or a cell array of strings
            %       threshMethod    (optional) 'or', 'and', determines thresholding logic
            %                        - Default: 'or'
            %
            %   Ouputs
            %   
            %       data            Updated field with .tf field for each color specified, which 
            %                       has the 'transfected' or gated + cells given the threshold
            %                           (note - if no threshold is supplied, data is unchanged.)
            %       stats           A struct containing population stats in the following fields,
            %                       each which is an H x W x numel(channels) matrix (where H and W
            %                       are the height and width of the data struct). 
            %                           all10s          10th percentile
            %                           allMedians      Median (50th percentile)
            %                           all90s          90th percentile
            %                           allMeans        Geometric mean 
            %                           allSDs          Geometric standard deviation
            %                           ratioMedians    (not implemented)
            %                           ratioMeans      (not implemented)
            %                           ratioSDs        (not implemented)
            %       
            % Written by Ross Jones
            %   Weiss Lab, MIT
            %   2016-05-27
            %   
            % Updates:
            %   2016-12-06: Added thresholding by average percentage above threshold (threshPct)
            
            % Check inputs
            [data, channels, dataType] = checkInputs(inputs);
            dataType = dataType{:}; % Not implemented multiple yet
            
            [H, W, D] = size(data);
            C = numel(channels);
            
            % If a single string is given for thresh channels, convert to a cell
            if (isfield(inputs, 'threshChan') && ischar(inputs.threshChan))
                inputs.threshChan = {inputs.threshChan};
            end
            
            if (~isfield(inputs, 'threshMethod'))
                inputs.threshMethod = 'or';
            end
            
            % Check if percentage-based thresholding and find percent if applicable
            if (isfield(inputs, 'threshPct') && logical(inputs.threshPct))
                % Check other required inputs are present
                assert(isfield(inputs, 'threshChan'), 'Must have ''threshChan'' input!')
                assert(isfield(inputs, 'thresh'), 'Must have ''thresh'' input!')
                
                % Find percentThresh
                percentPos = zeros(numel(data), 1);
                for d = 1:numel(data)
                    activeSet = threshold(data(d), inputs.threshChan, inputs.thresh, dataType, inputs.threshMethod);
                    percentPos(d) = sum(activeSet) / numel(activeSet) * 100;
                end
                percentThresh = mean(percentPos(:));
            else
                percentThresh = false;
            end
            
            all10s      =   zeros(H, W, C);
            allMedians  =   zeros(H, W, C);
            all90s      =   zeros(H, W, C);
            allMeans    =   zeros(H, W, C);
            allSDs      =   zeros(H, W, C);
            for i = 1:H
                for j = 1:W
                    
                    % Setup data vectors
                    chans = cell(1, numel(channels));
                    
                    % Collect data from all replicates
                    for k = 1:D
                        
                        % Threshold (ie gate transfected cells) if applicable
                        if all(isfield(inputs, {'threshChan', 'thresh'}))
                            if (percentThresh)
%                                 activeSet = false(size(data(i, j, k).(channels{1}).(dataType))); 
                                
                                % Add together fluorescence^2 from all channels and find the pct
                                % highest total expressing cells
                                % - We do this because it essentially takes all cells above some 
                                %   circular/spherical/etc limit surrounding the non-TF cells.
                                %   (Plot X^2 + Y^2 + Z^2 + ... > TH^2)
                                %    - Taking the top P percent from all channels will give many 
                                %      non-transfected cells
                                %    - Using multiplication to find ideal output would exclude 
                                %      single-pos cells (Plot X*Y > TH --> Y > TH / X)
                                %    - Using addition would bias towards single-pos cells
                                %      (Plot X + Y > TH --> Y > TH - X)
                                %    - Using the highest val in any given channel takes a square/
                                %      cube/etc from around cells (Plot Y > 1 & X > 1)
                                fluorScore = zeros(size(data(i, j, k).(channels{1}).(dataType)));
                                for ch = 1:numel(inputs.threshChan)
                                    chanData = data(i, j, k).(inputs.threshChan{ch}).(dataType);
                                    fluorScore = fluorScore + chanData.^2 .* sign(chanData);
                                end
                                pctThreshVal = prctile(fluorScore, 100 - percentThresh);
                                activeSet = (fluorScore > pctThreshVal);
                            else
                                activeSet = threshold(data(i, j, k), inputs.threshChan, inputs.thresh, dataType, inputs.threshMethod);
                            end
                        end
                        
                        % Iterate over replicates to collate data
                        for ch = 1:numel(channels)
                            
                            % Pull out data
                            chanData = data(i, j, k).(channels{ch}).(dataType);
                            
                            % If applicable, apply thresholding and save thresholded population
                            if all(isfield(inputs, {'threshChan', 'thresh'}))
                                chanData = chanData(activeSet);
                                data(i, j, k).(channels{ch}).tf = chanData;
                            end
                            
                            % Collect replicates - stats are calculated on combined data
                            chans{ch} = [chans{ch}; chanData];
                        end
                    end
                    
                    
                    % Calculate stats for each channel
                    for ch = 1:numel(channels)
                        
                        % Calculate geometric means, medians, and standard devs
                        pct = prctile(chans{ch}, [10, 50, 90]);
                        gm = geomean(chans{ch}(chans{ch} > 0));
                        sd = geostd(chans{ch}(chans{ch} > 0));
                        
                        all10s(i, j, ch) = pct(1);
                        allMedians(i, j, ch) = pct(2);
                        all90s(i, j, ch) = pct(3);
                        allMeans(i, j, ch) = gm;
                        allSDs(i, j, ch) = sd;
                    end
                end
            end
            
            % Check requested output stats
            outStats = {'all10s', 'allMedians', 'all90s', 'allMeans', 'allSDs'};
            if isfield(inputs, 'statTypes')
                
                % Convert single string entries to a cell
                if ischar(inputs.statTypes)
                    inputs.statTypes = {inputs.statTypes};
                end
                
                % Wean list
                outStats = intersect(inputs.statTypes, outStats);
            end
            
            % Create output stats struct
            stats = struct();
            for st = 1:numel(outStats)
                stats.(outStats{st}) = eval(outStats{st});
            end
            
%             % Get ratio of means to mimic
%             ratioMedians = zeros(size(allMedians));
%             ratioMeans = zeros(size(allMeans));
%             ratioSDs = zeros(size(ratioMeans));
%             for i = 1:H
%                 for j = 1:W
%                     for k = 1:numel(channels)
% 
%                         % Control (Gal4 only) data is in first index
%                         ratioMeans(i, j, k) = allMeans(i, j, k) ./ allMeans(1, 1, k);
%                         ratioMedians(i, j, k) = allMedians(i, j, k) ./ allMedians(1, 1, k);
% 
%                         % New SD must be calculated w/ formula s = e^sqrt(log(s1)^2 + log(s2)^2) since geostd
%                         ratioSDs(i, j, k) = exp(sqrt(log(allSDs(i, j, k)).^2 + log(allSDs(1, 1, k)).^2));
%                     end
%                 end
%             end
            
            
            % --- Helper Function --- %
            
            
            function [data, channels, dataType] = checkInputs(inputs)
                
                % Ensure necessary inputs are present
                reqFields = {'data', 'channels', 'dataType'};
                hasFields = isfield(inputs, reqFields);
                assert(all(hasFields), 'Missing field: %s\n', reqFields{~hasFields})
                
                % Extract data
                data = inputs.data;
                channels = inputs.channels;
                dataType = inputs.dataType;
                
                % Checks the inputs to make sure they are valid
                validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                if ischar(channels)
                    % Convert to cell if a single string
                    channels = {channels};
                end
                for n = 1:numel(channels)
                    validatestring(channels{n}, fieldnames(data(1)), mfilename, 'channel', 2);
                end
                
                % Convert dataType to cell if just a string
                if ~iscell(dataType)
                    dataType = {dataType};
                end
                
                % Check if more than one dataType is given
                if numel(dataType) > 1
                    assert(numel(dataType) == numel(channels), ...
                        'If more than one dataType given, must match number of channels')
                end
                
                % Validate dataTypes
                for n = 1:numel(dataType)
                    validatestring(dataType{n}, fieldnames(data(1).(channels{1})), mfilename, 'dataType', 3);
                end
            end
            
            
            function activeSet = threshold(threshData, threshChan, thresh, dataType, threshMethod)
                % Thresholds the given data on the given channel(s) with the
                % given threshold(s) on the given dataType.
                %
                %   Optional: threshMethod - 'or' (default), 'and' 
                %       Determines the logic of counting thresholds
                %
                % Returns a logical vector indicating cells that pass the thresh

                if ~exist('threshMethod', 'var')
                    threshMethod = 'or';
                end
                
                % Pre-allocate activeSet as a full-sized FALSE vector for OR
                % logic and a TRUE vector for AND logic
                switch threshMethod
                    case 'and'
                        activeSet = true(size(threshData.(threshChan{1}).(dataType)));
                    otherwise
                        activeSet = false(size(threshData.(threshChan{1}).(dataType)));
                end
                
                % For each thresholding channel, add passing cells to the activeSet
                for tc = 1:numel(threshChan)
                    if (numel(thresh) > 1)
                        assert(numel(thresh) == numel(threshChan), ...
                            'If more than one thresh, must equal # of thresholded channels')
                        currThresh = thresh(tc);
                    else
                        currThresh = thresh;
                    end
                    switch threshMethod 
                        case 'and'
                            activeSet = (activeSet & (threshData.(threshChan{tc}).(dataType) > currThresh));
                        otherwise
                            activeSet = (activeSet | (threshData.(threshChan{tc}).(dataType) > currThresh));
                    end
                end
            end
            
        end
        
        
        function data = mixtureModel(data, channels, dataType, numPop, par, xrange)
            % Fits the given data with a gaussian mixture model
            %
            %   Parameters:
            %       data        Struct: Standard data structure
            %       channels    Cell: A list of names of channels to fit on
            %       dataType    String: The type of data to fit (default = 'raw')
            %       numPop      Int: The number of populations to fit (default = 2)
            %       par         Logical: If true, uses parallel computing (default = false) 
            %       xrange      The range to use for pdf display (default = linspace(-2, 5, 1e5))
            %                       Note: Data is logicle transformed, then pdf is on log scale
            %
            %   Ouputs:
            %       data        Updated data struct. Each fitted channel now has sample means 
            %                   (.mus) and probability denisty functions (.pdfs) as fields.
            
            % Check existence of optional inputs, assign defaults
            if (~exist('dataType', 'var'))
                dataType = 'raw';
            end
            if (~exist('numPop', 'var'))
                numPop = 2;
            end
            if (~exist('par', 'var'))
                par = false;
            end
            if (~exist('xrange', 'var'))
                xrange = linspace(-2, 5, 1e5);
            end
            
            % Check inputs
            checkInputs(data, channels, dataType, numPop, par);
            numPop = round(numPop);
            
            % Find data sizes
            [H, W, D] = size(data);
            C = length(channels);
            
            
            % For tracking progress
            totalProgress = H * W * D * C;
            currProgress = 0;
            percentile = 10;
            
            % Set fitting parameters
            options = statset('MaxIter', 500);
            parameters = {
                'Start', 'plus', ...            % Use k-means++ to seed initial guesses
                'Replicates', 10, ...           % 10 replicates to avoid local minima
                'CovarianceType', 'full', ...   % Covariance of each component is not independent
                'SharedCovariance', false};     % Populations do not share covariance
            if (par)
                % In order to maximize parallel processing, ensure data struct is linear
                parfor i = 1:numel(data)
                    for channel = channels
                        
                        % Log transform the data first
                        % --> perhaps better to do hyperlog or biexp tranform?
                        d = data(i).(channel{:}).(dataType);
%                         logData = log10(d(d > 1e-2));
%                         badData = (isinf(logData) | isnan(logData));
%                         logData = logData(~badData);
                        logData = Transforms.lin2logicle(d);
                        
%                         xpdf = sort(logData(1:10:end));
%                         [~, pval] = hartigansDipSignifTest(xpdf, 5000);
%                         if (pval > 0.05)
%                             warning('No significant multi-modality detected, p = %.2f\n', pval)
%                         end
                        
                        % Fit Gaussian mixture model
                        gmm = fitgmdist(logData, numPop, 'Options', options, parameters{:}); %#ok<PFBNS>
                        
                        % Extract population means and probability density function (log space) 
                        data(i).(channel{:}).mus = sort(gmm.mu);
                        data(i).(channel{:}).pdfs = gmm.pdf(reshape(xrange, [], 1));
                    end
                    fprintf(1, 'Worker %d has finished\n', i);
                end
            else
                for i = 1:numel(data)
                    for channel = channels

                        % Log transform the data first
                        % --> perhaps better to do hyperlog or biexp tranform?
                        d = data(i).(channel{:}).(dataType);
                                logData = log10(d(d > 1e-2));
                                badData = (isinf(logData) | isnan(logData));
                                logData = logData(~badData);
%                         logData = Transforms.lin2logicle(d);

%                         xpdf = sort(logData(1:10:end));
%                         [~, pval] = hartigansDipSignifTest(xpdf, 5000);
%                         if (pval > 0.05)
%                             warning('No significant multi-modality detected, p = %.2f\n', pval)
%                         end
                        
                        % Fit Gaussian mixture model
                        gmm = fitgmdist(logData, numPop, 'Options', options, parameters{:});

                        % Extract population means and probability density function (log space) 
                        data(i).(channel{:}).mus = sort(gmm.mu);
                        data(i).(channel{:}).pdfs = gmm.pdf(reshape(xrange, [], 1));

                        % Update progress
                        currProgress = currProgress + 1;
                        pctProgress = currProgress / totalProgress * 100;
                        if (pctProgress > percentile)
                            fprintf(1, '%d%% Complete\n', percentile);
                            percentile = percentile + 10;
                        end
                    end
                end
            end
            
            
            % --- Helper Function --- %
            
            
            function checkInputs(data, channels, dataType, numPop, par)
                % Checks the inputs to make sure they are valid
                
                validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                validateattributes(channels, {'cell', 'char'}, {}, mfilename, 'channels', 2);
                if ischar(channels), channels = {channels}; end
                if (size(channels, 1) > 1 && size(channels, 2) > 1)
                    error('channels cannot be a 2D array of channel names');
                end
                if (~isequal(intersect(fieldnames(data(1)), channels), reshape(sort(channels), [], 1)))
                    error('Given channels not found in data');
                end
                validatestring(dataType, fieldnames(data(1).(channels{1})), mfilename, 'dataType', 3);
                validateattributes(numPop, {'numeric'}, {'positive'}, mfilename, 'numPop', 4);
                validateattributes(par, {'numeric', 'logical'}, {}, mfilename, 'par', 5);
            end
            
        end
        
        
        function [x1, x2, p1, p2, yfit1, yfit2, rsq_both] = LSB_fit(constitutive, reporter, doLog, skip)
            % Written by Jeremy Gam, Sept 4th 2014
            % This script is for analyzing single miRNA low sensor repression from
            % transfection data.
            %
            % Take in vectors for transfection marker expression (x), and sensor
            % expression (y)
            %
            % We assume two regimes in the data, one where the miRNA
            % target is repressed and increases in transfection marker expression
            % result in low increases in sensor expression and one where boeth marker
            % and sensor increase roughly 1:1. First we rank order data based on
            % transfection marker expression, log transform, split the marker
            % expression vector into all possible halves, and calculate SSresid for
            % linear fits for both halves and add them together. Take the split with
            % lowest total SSresid and you then have the piecewise fit and can find
            % points like the intercept.
            %
            %   Modified 2015-04-02 by Ross Jones
            %       - Constrained fits to only positive slopes
            %       - Added options for converting to log and passing skip as an arg
            %       - Added some input checking
            %       - Handling row vectors
            %       - Some readability improvements
            
            
            % --- Check Inputs --- %
            
            % Must have > 10 points to fit to
            if (length(constitutive) < 10 || length(reporter) < 10)
                error('Length of input vectors must be > 10')
            end
            
            % Check if doLog was passed, if not, default it to false
            if (~exist('doLog', 'var'))
                doLog = false;
            end
            
            % If dictated, the first step is to transform data to log space
            %   The distribution is more 'normal' in log than linear
            if (doLog)
                constitutive = log10(constitutive);
                reporter = log10(reporter);
            end
            
            % Check if skip was passed 
            if (~exist('skip', 'var'))
                % For speed, you can only calculate SSresid for every "skip" cells. 
                %   For skip = 100, only 1% of cells are used
                skip = 100;
            end
            
            % --- Sort Data --- %
            
            % Build a data matrix
            %   The first step ensures the data are in column vectors
            %   Then they are appended next to each other in an Nx2 matrix
            constitutive = reshape(constitutive, [], 1);
            reporter = reshape(reporter, [], 1);
            data = [constitutive, reporter];

            % Sort the data by transfection marker
            data = sortrows(data, 1);
%             display(data)
%             plot(data(:,1),data(:,2),'o')
            
            % --- Do Fitting --- %
            
            %{
            This is the linear regression example from MathWorks

            load count.dat
            x = count(:,1);
            y = count(:,2);

            p = polyfit(x,y,1);
            display(p)

            yfit = polyval(p,x);
            yresid = y - yfit;
            SSresid = sum(yresid.^2);
            SStotal = (length(y)-1) * var(y);
            rsq = 1-SSresid/SStotal;
            display(rsq)
            %}
            
            % Split cells into two groups by low/high flourescence of marker into every
            % possible set of halves and iteratively calc SSresid - sum of squares of
            % residuals (SSresid is metric for how good the fit is). Choose the ideal
            % set of two groups as lowest SSresid.
            N = length(constitutive);
            update = 5;
            progress = update;
            SSresid = nan(N, 1);
            for i = 2 : skip : N - 1 
                % Loop through every possible division into low/high
                % see note above about skip
                
                first_half = data(1:i, :);
                second_half = data(i:end, :);
%                 display(first_half)
%                 display(second_half)
                
                % Fit first half
                x1 = first_half(:, 1); % marker flourescence for low marker half
                y1 = first_half(:, 2); % reporter flourescence for low marker half
                p1 = polyfit(x1, y1, 1);
                yfit1 = polyval(p1, x1);
                yresid1 = y1 - yfit1;
                if (doLog)
                    % Non-linear operation, so need to transform back to linear scale
                    yresid1 = 10.^(yresid1);
                end
                SSresid1 = sum(yresid1.^2);
                
                % Fit second half
                x2 = second_half(:, 1); % marker flourescence for high marker half
                y2 = second_half(:, 2); % reporter flourescence for high marker half
                p2 = polyfit(x2, y2, 1);
                yfit2 = polyval(p2, x2);
                yresid2 = y2 - yfit2;
                if (doLog)
                    % Non-linear operation, so need to transform back to linear scale
                    yresid2 = 10.^(yresid2);
                end
                SSresid2 = sum(yresid2.^2);
                
                % Add residuals together and store
                SSresid(i) = SSresid1 + SSresid2;
                
                % Check if either line has a negative slope - we want pos slopes
                if (p1(1) < 0 || p2(1) < 0)
                    SSresid(i) = inf;
                end
                
                % Display progress
                %   Only if skip is set high, meaning we expect a lot of data
                if (skip >= 100 && i / N * 100 > progress)
                    fprintf(1, '%d%% complete\n', progress);
                    progress = progress + update;
                end
            end
            
            % Make sure we don't take the boundaries as the ideal division point
            SSresid([1, end]) = inf; 

            % Find the minimum of SSresid and recalculate the fit
%             display(SSresid)
            [~, idx] = min(SSresid);
            
            % Assign optimal halves
            first_half = data(1:idx, :);
            second_half = data(idx:end, :);
%             display(first_half)
%             display(second_half)
            
            % Refit first half
            x1 = first_half(:, 1);
            y1 = first_half(:, 2);
            p1 = polyfit(x1, y1, 1);
            yfit1 = polyval(p1, x1);
            yresid1 = y1 - yfit1;
            %SSresid1 = sum(yresid1.^2);
            %SStotal1 = (length(y1)-1) * var(y1);
            %rsq_adj1 = 1 - SSresid1/SStotal1 * (length(y1)-1)/(length(y1)-length(p1));
            
            % Refit second half
            x2 = second_half(:, 1);
            y2 = second_half(:, 2);
            p2 = polyfit(x2, y2, 1);
            yfit2 = polyval(p2, x2);
            yresid2 = y2 - yfit2;
            %SSresid2 = sum(yresid2.^2);
            %SStotal2 = (length(y2)-1) * var(y2);
            %rsq_adj2 = 1 - SSresid2/SStotal2 * (length(y2)-1)/(length(y2)-length(p2));
            
%             figure()
%             plot(constitutive, reporter, x1, yfit1, x2, yfit2)
            
            % Sum the square of the residuals
            SSresid_both = sum(yresid1.^2) + sum(yresid2.^2);
            SStotal_both = (length(y1) - 1) * var(y1) + (length(y2) - 1) * var(y2);
            rsq_both = 1 - SSresid_both / SStotal_both;

            % Find and return metrics of interest
%             out = [p1, p2, rsq_both];
            
            % Transform fit to linear space if data was log transformed
            if (doLog)
                x1 = 10.^(x1);
                x2 = 10.^(x2);
                p1 = 10.^(p1);
                p2 = 10.^(p2);
                yfit1 = 10.^(yfit1);
                yfit2 = 10.^(yfit2);
                % rsq is already made linear earlier
            end
        end

        
        function data = ronBin(data, constitutive, outputs)
            % Creates Ron's style of bins for the given data, adding the 'bins' 
            % and 'means' fields to the data struct. 
            %
            % Inputs
            %   'data'          A data struct like that generated by the openFiles() method.
            %   'constitutive'  A string name of the constitutive reporter channel. Must match 
            %                   a field in the data struct with the subfield 'raw'.
            %   'outputs'       Similar to constitutive, but for output markers. Can be a cell 
            %                   array of channel names or a single channel name.
            %
            % Outputs
            %   'means'         A substruct with fields which are the names of the the 
            %                   output channels given as input. Each of these fields
            %                   points to a 1x1000 array of doubles containing the means
            %                   of each bin. 
            %   'edges'         Similarly, a substruct with fiels for the output reporters.
            %                   The bin edges are given as a double array, and correspond 
            %                   with the data in 'means'.
            
            % Check inputs
            validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
            validateattributes(constitutive, {'char'}, {}, mfilename, 'constitutive', 2);
            validateattributes(outputs, {'char', 'cell'}, {}, mfilename, 'outputs', 3);
            if (~iscell(outputs))
                outputs = {outputs};    % This just makes things less messy throughout
            end
            numOutputs = numel(outputs);
            
            % Define bins
            AXES_MIN = 1;
            AXES_MAX = 2.5e+5;
            NUM_BINS = 100;
            binEdges = logspace(log10(AXES_MIN), log10(AXES_MAX), NUM_BINS + 1);
            
            % Iterate through and bin the data
            for i = 1:length(data)
                sampleConst = data(i).(constitutive).raw;                
                
                % Find what bins each data point is in by looking through the 
                % constitutive reporter.
                [~, whichBin] = histc(sampleConst, binEdges);
                
                % Find means of output bins
                for j = 1:numOutputs
                    sampleOut = data(i).(outputs{j}).raw;
                    binMean = zeros(1, NUM_BINS);
                    for k = 2:NUM_BINS
                        flagBinMembers = (whichBin == k);
                        binMembers     = sampleOut(flagBinMembers);
                        binMean(k)     = mean(binMembers);
                    end
                                        
                    % Select out those that are greater than zero
                    flagBinMeans = (binMean > 0);
                    nonzeroBinMeans = binMean(flagBinMeans);
                    nzBE = binEdges(1:NUM_BINS);
                    nonzeroBinEdges = nzBE(flagBinMeans);
                    
                    % Smooth the bins
                    %yy = smooth(realBinEdges, realBinMeans, 0.1,'rloess');
                    smoothedMeans = smooth(nonzeroBinEdges, nonzeroBinMeans);
                    
                    data(i).means.(outputs{j}) = smoothedMeans; 
                    data(i).edges.(outputs{j}) = nonzeroBinEdges;
                end
            end
            
            fprintf(1, 'Finished binning data\n');
        end
    end
    
end