classdef Transforms < handle
	% Compilation of data transforms from the Weiss Lab flow cytometry MATLAB repository
	%
	%	Methods are implemented as static functions, so they can be called directly
	%	or after creating a handle to this class.
	%
	%		Ex:
	%		tfData = Transforms.lin2logicle(data);
	%
	%		Ex2:
	%		tf = Transforms();
	%		tfData = tf.logicle2lin(data);
	%
    %   Functions:
    %
    %       out = lin2logicle(S)
    %       out = logicle2lin(X)
	%		beadVals = getBeadVals(beadType, beadLot, channels)
    %       data = fcs2MEF(data, channelFits, dataType)
    %       [channelFits, beadDir] = calibrateMEF(beadsFilename, channels, showPlots)
    %       data = fcs2ABC(data, channelFits, dataType)
	%		channelFits = calibrateABCfromMedians(beadsFilenames, channel, hasBlank)
    %       [channelFits, locs, peaks] = calibrateABC(beadsFilename, channel, largePeaks, locs, peaks)
    %
    % Written/Compiled by Ross Jones
    % Weiss Lab, MIT
    % Last updated 2017-08-10
	
	properties (Constant)
		
		MEF_CONVERSION_FACTOR = 2e3;
		
		CHANNEL_MAP = struct( ...
			'Cascade_Blue_A',	'MEVSB', ...
			'Pacific_Blue_A',	'MEBFP', ...
			'FITC_A',			'MEFL', ...
			'PE_A',				'MEPE', ...
			'PE_Texas_Red_A',	'MEPTR', ...
			'PE_TxRed_YG_A',	'MEPTR', ...	% Alt name for Koch LSRII-HTS2
			'PE_Cy5_5_A',		'MECY', ...
			'PE_Cy7_A',			'MEPCY7', ...
			'APC_A',			'MEAP', ...
			'APC_Cy7_A',		'MEAPCY7');
		
		BEAD_TYPES = {'RCP-30-5A'};
		
		RCP305A_LOTS_001 = {'AD04', 'AE01', 'AF01', 'AF02', ...
							'AH01', 'AH02', 'AJ01'};
		RCP305A_LOTS_002 = {'AA01', 'AA02', 'AA03', 'AA04', ...
							'AB01', 'AB02', 'AC01', 'GAA01-R'};
		
		RCP305A_VALS_001 = struct( ... % from http://www.spherotech.com/RCP-30-5a%20%20rev%20H%20ML%20071712.xls
			'MECSB',	[216, 464, 1232, 2940, 7669, 19812, 35474], ...
			'MEBFP',	[861, 1997, 5776, 15233, 45389, 152562, 396759], ...
			'MEFL',		[792, 2079, 6588, 16471, 47497, 137049, 271647], ...
			'MEPE',		[531, 1504, 4819, 12506, 36159, 109588, 250892], ...
			'MEPTR',	[233, 669, 2179, 5929, 18219, 63944, 188785], ...
			'MECY',		[1614, 4035, 12025, 31896, 95682, 353225, 1077421], ...
			'MEPCY7',	[14916, 42336, 153840, 494263], ...
			'MEAP',		[373, 1079, 3633, 9896, 28189, 79831, 151008], ...
			'MEAPCY7',	[2864, 7644, 19081, 37258]);
		
		RCP305A_VALS_002 = struct( ... % from http://www.spherotech.com/RCP-30-5a%20%20rev%20G.2.xls
			'MECSB',	[179, 400, 993, 3203, 6083, 17777, 36331], ...
			'MEBFP',	[700, 1705, 4262, 17546, 35669, 133387, 412089], ...
			'MEFL',		[692, 2192, 6028, 17493, 35674, 126907, 290983], ...
			'MEPE',		[505, 1777, 4974, 13118, 26757, 94930, 250470], ...
			'MEPTR',	[207, 750, 2198, 6063, 12887, 51686, 170219], ...
			'MECY',		[1437, 4693, 12901, 36837, 76621, 261671, 1069858], ...
			'MEPCY7',	[32907, 107787, 503797], ...
			'MEAP',		[587, 2433, 6720, 17962, 30866, 51704, 146080], ...
			'MEAPCY7',	[718, 1920, 5133, 9324, 14210, 26735]);
	end
	
	methods (Static)
		
		function out = lin2logicle(S)
			% lin2logicle(S) converts a vector of values from a linear scale to a logicle scale 
			%   out = lin2logicle(S) returns a vector of values corresponding to the values 
			%   in X transformed from a linear scale to a logicle one using
			%   the conversion function described by Parks, et al. "A New Logicle
			%   Display Method Avoids Deceptive Effects of Logarithmic Scaling for Low
			%   Signals and Compensated Data" Equation (5), where
			%     S - a vector of linear 'raw' values
			%     out - a vector of logicle converted values from S
			%
			%   Since Equation 5 cannot be explicitly solved for X, a spline is fitted to
			%   X and S data. 
			%
			%   Example:
			%       out = Transforms.lin2logicle(linspace(0,1000));
			%
			%   Written by
			%   Breanna Stillo
			%   bstillo@mit.edu
			%   Last Updated: 2014-10-14
			
			X = linspace(0, 4.5, 1000); % 4.5 = lin2logicle(2^18)
			Y = Transforms.logicle2lin(X);
			p = spline(Y, X);

			out = ppval(p, S);

		end
		
		
		function out = logicle2lin(X)
			% logicle2lin(X) converts a vector of values from a logicle scale to a linear scale 
			%   out = logicle2lin(X) returns a vector of values corresponding to the values 
			%   in X transformed from a logicle scale to a linear one using
			%   the conversion function described by Parks, et al. "A New Logicle
			%   Display Method Avoids Deceptive Effects of Logarithmic Scaling for Low
			%   Signals and Compensated Data" Equation (5), where
			%     X - a vector of logicle values
			%     out - a vector of linear converted values from X
			%
			%   Example:
			%       out = Transforms.logicle2lin(linspace(0,4));
			%
			%   Written by
			%   Breanna Stillo
			%   bstillo@mit.edu
			%   Last Updated: 2014-10-14
			%
			%   Edit 2015-02-09 (Ross) - Added comments and intuitive variable names
			%                            Major speed up by using algebraic conversion 
            %                            instead of symbolic toolbox
            %   Edit 2016-03-28 (Ross) - Made display breadth smaller w/ Bre's adjustment of 
            %                            DATA_MIN, which increases the relative size of the area
            %                            between 10^-2 and 10^2

			% Set values
			DATA_MAX = 2^18;         % (T) The maximum measurable data value ( = 26214 )
			DISPLAY_BREADTH = 4.5;   % (M) Breadth of the display in decades ( = lin2logicle(2^18) )
			DATA_MIN = -150;         % (r) Negative range reference value    ( lin2logicle(0) = -171 )

			% This gives a range for linearization around zero (W).
			linRange = (DISPLAY_BREADTH - log10(DATA_MAX / abs(DATA_MIN))) / 2;

			% p and W (linRange) are considered one parameter, p is introduced by the 
			% authors for compactness. Here we find p that solves their equivalence:
			%   w = 2p ln(p)/(p + 1)
			% Solved by WolframAlpha:
			P = linRange / (2 * lambertw(0.5 * exp(-linRange / 2) * linRange));

			% Find where the given data is out of the linear range.
			logX = (X >= linRange);

			% Compute and return the final linearized vector 
			out = toLin(X) .* logX - toLin(2 * linRange - X) .* (1 - logX);
			
			
			% --- Helper Functions --- %
			
			
			function converted = toLin(convert)
				% Does a final conversion
				
				converted = DATA_MAX .* 10^-(DISPLAY_BREADTH - linRange) .* ...
					(10.^(convert - linRange) - P^2 .* 10.^(-(convert - linRange) ./ P) + P^2 - 1);
			end
		end
		
        
		function out = lin2logicleMEF(S)
			% lin2logicleMEF(S) converts a vector of MEF values from a linear scale to a logicle scale 
			%   out = lin2logicleMEF(S) returns a vector of values corresponding to the values 
			%   in S transformed from a linear scale to a logicle one using
			%   the conversion function described by Parks, et al. "A New Logicle
			%   Display Method Avoids Deceptive Effects of Logarithmic Scaling for Low
			%   Signals and Compensated Data" Equation (5), where
			%     S - a vector of linear 'raw' values
			%     out - a vector of logicle converted values from S
			%
			%   Since Equation 5 cannot be explicitly solved for X, a spline is fitted to
			%   X and S data. 
			%
			%   Example:
			%       out = Transforms.lin2logicleMEF(linspace(0,1000));
			%
			% Written by
			% Ross Jones
			% jonesr18@mit.edu
			%
			%	Updated Log: 
			%		2017-09-06	Converted this from Transforms.lin2logicle()
			%					for use with MEF-converted data. 
			
			% Scale down MEF data to be comparable to raw data
			%	The yellow/green FPs are very bright and MEFL units will be ~2
			%	decades higher than fluorescence.
			S = S / Transforms.MEF_CONVERSION_FACTOR; 
			
			X = linspace(0, 4.5, 1000); % 4.5 = lin2logicle(2^18)
			Y = Transforms.logicle2lin(X);
			p = spline(Y, X);
			
			out = ppval(p, S);

		end
		
		
		function out = logicle2linMEF(X)
			% logicle2linMEF(X) converts a vector of values from a logicle scale to a linear scale 
			%   OUT = logicle2linMEF(X) returns a vector of values corresponding to the values 
			%   in X transformed from a logicle scale to a linear one using
			%   the conversion function described by Parks, et al. "A New Logicle
			%   Display Method Avoids Deceptive Effects of Logarithmic Scaling for Low
			%   Signals and Compensated Data" Equation (5), where
			%     X - a vector of logicle values
			%     out - a vector of linear converted values from X
			%
			%   Example:
			%       out = Transforms.logicle2linMEF(linspace(0,4));
			%
			% Written by
			% Ross Jones
			% jonesr18@mit.edu
			%
			%	Updated Log: 
			%		2017-09-06	Converted this from Transforms.logicle2lin()
			%					for use with MEF-converted data. 
			
			% Set values
			DATA_MAX = 2^18;			% (T) The maximum measurable data value ( ~2.6E+5 )
			DISPLAY_BREADTH = 4.5;		% (M) Breadth of the display in decades ( = lin2logicle(2^18) )
			DATA_MIN = -1.5e2;			% (r) Negative range reference value    ( logicle2lin(0) = -1.5E+2 )

			% This gives a range for linearization around zero (W).
			linRange = (DISPLAY_BREADTH - log10(DATA_MAX / abs(DATA_MIN))) / 2;

			% p and W (linRange) are considered one parameter, p is introduced by the 
			% authors for compactness. Here we find p that solves their equivalence:
			%   w = 2p ln(p)/(p + 1)
			% Solved by WolframAlpha:
			P = linRange / (2 * lambertw(0.5 * exp(-linRange / 2) * linRange));

			% Find where the given data is out of the linear range.
			logX = (X >= linRange);

			% Compute and return the final linearized vector 
			%	Note that the output is scaled back up to MEF-range units 
			out = Transforms.MEF_CONVERSION_FACTOR * toLin(X) .* logX - toLin(2 * linRange - X) .* (1 - logX);
			
			
			% --- Helper Functions --- %
			
			
			function converted = toLin(convert)
				% Does a final conversion
				
				converted = DATA_MAX .* 10^-(DISPLAY_BREADTH - linRange) .* ...
					(10.^(convert - linRange) - P^2 .* 10.^(-(convert - linRange) ./ P) + P^2 - 1);
			end
		end
		
		
		function beadVals = getBeadVals(beadType, beadLot, channels)
			% Returns the bead values for the given beads of the given lot in the given channels. 
			% 
			%	beadVals = getBeadVals(beadType, beadLot, channels)
			%
			%	Inputs 
			%		beadType		<char> The name of the type of bead (eg 'RCP-30-5A')
			%		beadLot			<char> The bead production lot (eg 'AH01')
			%		channels		<cell, char> The channels to extract bead
			%						values for (eg 'Pacific_Blue_A')
			%
			%	Outputs
			%		beadVals		<struct> A struct where bead unit names (eg 'MEF'  
			%						are the keys and the bead peaks are the values
			
			% Check inputs
			validBeadLots = [Transforms.RCP305A_LOTS_001, Transforms.RCP305A_LOTS_002];
			validatestring(beadType, Transforms.BEAD_TYPES, mfilename, 'beadType', 1);
			validatestring(beadLot, validBeadLots, mfilename, 'beadLot', 2);
			validateattributes(channels, {'char', 'cell'}, {}, mfilename, 'channels', 3);
			badChannels = setdiff(channels, fieldnames(Transforms.CHANNEL_MAP));
			assert(isempty(badChannels), ...
				'Channel not valid: %s', badChannels{:});
			
			% Convert channel names to MEF units
			requestedUnits = cell(1, numel(channels));
			for i = 1:numel(channels), requestedUnits{i} = Transforms.CHANNEL_MAP.(channels{i}); end
			
			% Collect values into struct
			beadVals = struct();
			switch beadType
				case {'RCP-30-5A'}
					switch beadLot
						case Transforms.RCP305A_LOTS_001
							for u = requestedUnits
								beadVals.(u{:}) = Transforms.RCP305A_VALS_001.(u{:});
							end
						case Transforms.RCP305A_LOTS_002
							for u = requestedUnits
								beadVals.(u{:}) = Transforms.RCP305A_VALS_002.(u{:});
							end
					end
			end
		end

		
		function data = fcs2MEF(data, channelFits, dataType)
			% Converts fcs data to MEF units
			%	For conversion information, bead types/lots, valid channels, and
			%	MEF unit names, see Transforms.getBeadVals().
            %
            %   Inputs
            %
			%       data            A standard data struct. The channels must match those in channelFits. 
            %
            %       channelFits     (Output from calibrateMEF())
            %                       An Nx2 table of slopes/intercepts for linear conversion of 
            %                       N channels to MEF equivalents. Rows are labeled with 
            %                       channel names and columns are labeled with slope / zero
            %
            %       dataType        The population to convert to MEFs (eg 'scComp' or 'tf')
            %                        - must be in data!
			%
			%	Outputs
			%
			%		data			Data struct with MEF units added to all converted 
			%						channels in a new field called 'mef'.
			%						Also adds new gate 'nneg' - logical array flagging
			%						non-negative input data.
			%
			%   Written by
			%   Breanna Stillo
			%   bstillo@mit.edu
			%
            %   Update 2016-04-23 by Ross Jones
            %       Uses standard struct interface and uses channelFits directly rather than
            %       loading the saveFile
			
            % Process inputs
			channels = checkInputs_fcs2MEF();
                        
			for i = 1:numel(data) %#ok<ALIGN>
                
				if isempty(data(i).(channels{1})) % Some FlowData tcData will be empty 
					continue
				end
				
				% Identify non-negative cells
				nonNeg = true(size(data(i).(channels{1}).(dataType)));
				
				for ch = 1:numel(channels)
					
					% Extract channel data
                    chData = data(i).(channels{ch}).(dataType);
					chNeg = (chData < 0);
					
					% Identify non-neg data for gating
					nonNeg = (nonNeg & ~chNeg);
                    
					% Convert data - take abs value and flip neg data so it is not complex 
					mefData = 10.^polyval(channelFits.(channels{ch}), log10(chData));
					mefData(chNeg) = -abs(mefData(chNeg));
% 					mefData = polyval(channelFits.(channels{ch}), (chData));
					
                    % Do conversion
                    data(i).(channels{ch}).mef = mefData;
				end
				
				data(i).gates.nneg = nonNeg;
            end
			
			
			% --- Helper Functions --- %
			
            
			function channels = checkInputs_fcs2MEF()
				
				validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                validateattributes(channelFits, {'table'}, {}, mfilename, 'channelFits', 2);
                
                channels = channelFits.Properties.VariableNames;
                badChannels = setdiff(channels, fieldnames(data(1)));
                assert(isempty(badChannels), ...
                    'Channel not in data: %s\n', badChannels{:});
                
                validatestring(dataType, fieldnames(data(1).(channels{1})), mfilename, 'dataType', 3);
			end

        end
		
        
        function [channelFits, beadDir] = calibrateMEF(beads, channels, options)
            % calibrateMEF creates MEF fits for an .fcs file of a particular bead
            % sample chosen by the user. Note this is for rainbow calibration beads, not 
            % other antibody-tagged beads (see calibrateABC()). 
            %
            % A subfolder is created with the name of the bead file where associated calibration 
            % plots are stored and the calibration table corresponding to MEF fits.
            %
            %   Inputs: 
            %   
			%		beads			<struct> A struct with the following fields:
			%			filename		The .fcs file containing bead data
			%							corresponding with this experiment
			%			type			The name of the type of bead (eg 'RCP-30-5A')
			%			lot				The bead production lot (eg 'AH01')
			%		channels		<cell> Cell array of channel names to convert to MEF.
			%		options			<cell> Contains optional string flags:
			%			showPlots		Flag to show fitted plots
			%			nonLinear		Flag to *not* force linear bead fits
			%							--> This may be desired because theoretically
			%							the MEF-fluorescence relationship is linear,
			%							but is often slightly off. 
			%			fitGaussian		Flag to use multi-dimensional gaussian
			%							fitting rather than the single-dim 
			%							peak-finding algorithm
			%			peak_X			Forces peak X (where X is replaced with
			%							a number) to be the highest bead peak,
			%							rather than having the function iterate
			%							through all possible high peaks automaticallly
			%			noFilter		Do not filter out peaks by height
            %
            %   Outputs:
            %
            %       channelFits     <table> An Nx2 table of slopes/intercepts for  
            %                       linear conversion of N channels to MEF equivalents. 
            %                       Rows are labeled with channel names and columns are 
			%						labeled with slope / zero.
            %       beadDir         <char> The folder name for storing bead fits 
            %
            %   Written by
            %   Breanna Stillo & Ross Jones
            %   bstillo@mit.edu
			%	jonesr18@mit.edu
            %   Last Updated: 2017-09-16 by Ross Jones
			
            % Check inputs, initialize data struct
            [beadData, beadVals] = checkInputs_calMEF();
			MEF_units = fieldnames(beadVals)'; % Start in same order as channels
			
			% Setup new directory for fitting files/figs
			beadDir = ['Calibration', filesep];
			if ~exist(beadDir, 'file')
				mkdir(beadDir)
			end
			
			% Check if these exact beads have already been used for calibraiton
			% --> If so, load them and immediately return!
			beadSaveName = [beads.date, '_', beads.type, '_', beads.lot, '_', beads.cytometer];
			fitsFilename = [beadDir, beadSaveName '_MEF_Fits.mat'];
			if exist(fitsFilename, 'file')
				fprintf(1, 'Loading pre-computed bead fits\n');
				load(fitsFilename);
				return
			end
			
			% Find max number of peaks to look for
			% --> Some channels only have beads above a certain fluorescence
			BEAD_PEAKS = 8;
			MAX_PEAKS = inf;
			for chID = 1:numel(MEF_units)
				MAX_PEAKS = min(MAX_PEAKS, numel(beadVals.(MEF_units{chID}))); 
			end
			
			% Extract bead data and compute single-channel histograms
			[extractedBeadData, figFits, counts, centers] = extractBeadData(beadData, ismember('showPlots', options));
			
			% Find peaks using either Gaussian fits or MATLAB's findpeaks() function			
			if ismember('fitGaussian', options)
				warning('off', 'stats:gmdistribution:IllCondCov')
				warning('off', 'stats:gmdistribution:FailedToConvergeReps')
				means = fitGaussians(extractedBeadData);
			else
				means = fitFindPeaks(extractedBeadData);
			end
			
			% Filter means
			thresh = 8;
			filteredMeans = filterMeans(means, counts, centers, thresh);
			
			numPeaks = size(filteredMeans, 1) - 1;			
			res_min = inf;
			opt = optimoptions('lsqcurvefit', 'Display', 'off');
			for highestPeak = 1 : (MAX_PEAKS - numPeaks + 1)
				% Its a bit simpler to index peaks from the highest down,
				% since all beads have the largest populations, but some
				% don't have smaller ones. By flipping so that low indexes 
				% are high beads and vice-versa, we don't have to adjust 
				% the index for each bead population up or down. 
				
				highestPeak_true = BEAD_PEAKS - highestPeak + 1;
				fprintf(1, 'Trying highest peak: %d\n', highestPeak_true);
				
				% Iterate over each channel and add up residuals
				% Reset rolling parameters each time for cleanliness
				ress = zeros(1, numel(channels));
				fits = zeros(2, numel(channels));
				for chID = 1:numel(channels)

					% Extract vals 
					MEF = flipud(beadVals.(MEF_units{chID})');

					% Identify which bead pops to query
					beadPops = highestPeak : (highestPeak + numPeaks - 1);
					pointsMean = filteredMeans(2:end, chID); % Means already log10 transf
					pointsMEF = flipud(log10(MEF(beadPops)));
% 						pointsMEF = flipud((MEF(beadPops)));

					% Fit bead fluorescence to MEF values
					linearFunc = @(p, x) p(1) * x + p(2);
					fitNL = lsqcurvefit(linearFunc, [1; 1], pointsMean, pointsMEF, [], [], opt);
					if ismember('nonLinear', options)
						fit = fitNL;
					else
						fit = lsqcurvefit(@(p, x) linearFunc([1; p(2)], x), [1; 1], pointsMean, pointsMEF, [], [], opt);
					end
					fits(:, chID) = reshape(fit, [], 1);

					% Measure and collate residuals
					yResid = abs(pointsMEF - polyval(fitNL, pointsMean));
					ssResid = sum(yResid.^2);
					ssTotal = (length(pointsMEF) - 1) * var(pointsMEF);
					ress(chID) = ssResid / ssTotal;
				end
				res = min(ress);

				fprintf(1, 'Residuals: %.3f\n', res);
				
				% Determine if this is the best fit
				if ismember(sprintf('peak_%d', highestPeak_true), options)
					res = -inf;
					fprintf(1, 'Forcing Peak %d to be highest\n', highestPeak_true);
				end
				if (res < res_min)
					res_min = res;
					ress_min = ress;
					numPeaks_min = numPeaks;
					highestPeak_min = BEAD_PEAKS - highestPeak + 1;
					beadPops_min = BEAD_PEAKS - beadPops + 1;
					fits_min = fits;
					means_min = 10.^filteredMeans(2:end, :);
% 						gmmMeans_min = means(2:end, :);
				end
			end
			
			% Extract finalized value
			rsq			= 1 - ress_min;
			numPeaks	= numPeaks_min;
			highestPeak = highestPeak_min;
			beadPops	= beadPops_min;
			fits		= fits_min;
			meansLin	= means_min;
			meansBxp	= Transforms.lin2logicle(meansLin);
			
			fprintf(1, '\nFinished fitting! # Peaks = %d, Highest = Peak #%d\n', numPeaks, highestPeak);
			
			% Show final calibration plots if requested
			if ismember('showPlots', options)
				for chID = 1:numel(channels)
					% Find histogram bins that closest match GMM means
					locsRaw = zeros(1, size(means, 1));
					for muID = 1:size(means, 1)
						[~, locsRaw(muID)] = min(abs( ...
								centers.(channels{chID}) - means(muID, chID)));
					end
					locsFilt = zeros(1, numPeaks);
					for muID = 1:numPeaks
						[~, locsFilt(muID)] = min(abs( ...
								centers.(channels{chID}) - meansBxp(muID, chID)));
					end
					
					% Plot peaks on histogram
					figure(figFits.(channels{chID}));
					ax1 = subplot(2, 1, 1); hold(ax1, 'on');
					plot(ax1, means(:, chID), counts.(channels{chID})(locsRaw), '*b')
					plot(ax1, meansBxp(:, chID), counts.(channels{chID})(locsFilt), '*r')
					title(strrep(channels{chID}, '_', '\_'))
					ylabel('Count')
					xlabel('Fluorescence (AFU)')
					
					% Plot MEF fits, fit line, and squared sum of residuals
					MEF = fliplr(beadVals.(MEF_units{chID}));
					MEF = fliplr(MEF(BEAD_PEAKS - beadPops + 1));
					peakIntensity = meansLin(:, chID);
					
					% Plot fits
					ax2 = subplot(2, 1, 2);
					plot(ax2, ...
						peakIntensity, MEF, 'ob', ...
						peakIntensity, 10.^polyval(fits(:, chID), log10(peakIntensity)), 'r-')
					ax2.XScale = 'log';
					ax2.YScale = 'log';
					ax2.FontSize = 14;
					
					% Compute R^2 from residuals, display in title
					title(sprintf('Scale Factor = %.2f | R^2 = %.3f', ...
							10.^fits(2, chID), rsq(chID)), 'fontsize', 14)
					xlabel(ax2, 'Fluorescence')
					ylabel(ax2, MEF_units{chID});
					
					% Save figures to bead directory
					saveas(figFits.(channels{chID}), ...
								[beadDir, beadSaveName, '_', MEF_units{chID}, '_Fit.fig']); 
				end
			end
			
			% Convert fits to a table
			channelFits = array2table(fits);
			channelFits.Properties.RowNames = {'slope', 'zero'};
			channelFits.Properties.VariableNames = channels;
			
			% Save channel fits (and figures if applicable)
			save(fitsFilename, 'channelFits')
			
			
			% --- Helper Functions --- %
			
			
			function [beadData, beadVals] = checkInputs_calMEF()
				
				validateattributes(beads, {'struct'}, {}, mfilename, 'beadsFilename', 1);
				validFields = {'filename', 'type', 'lot', 'date', 'cytometer'};
				badFields  = setdiff(fieldnames(beads), validFields);
				assert(isempty(badFields), 'Field not valid: %s\n', badFields{:})
				
				beadsFilename = beads.filename;
				assert(logical(exist(beadsFilename, 'file')), ...
					'File does not exist: %s\n', beadsFilename)				
				
				% Type and lot are checked in Transforms.getBeadVals();
				
				validateattributes(channels, {'char', 'cell'}, {}, mfilename, 'channels', 4);
				
				% Convert single chars to cells for simplicity in processing
				if ~iscell(beadsFilename), beadsFilename = {beadsFilename}; end
				if ~iscell(channels), channels = {channels}; end
				
				% Import beads data
				beadData = FlowAnalysis.openFiles(beadsFilename{:});
				
				% Check channels are valid
				badChannels = setdiff(channels, fieldnames(beadData(1)));
				assert(isempty(badChannels), ...
					'Channel does not exist: %s\n', badChannels{:});
				
				% Get bead MEF values
				beadVals = Transforms.getBeadVals(beads.type, beads.lot, channels);
			end
			
			
			function [extractedBeadData, figFits, counts, centers] = extractBeadData(beadData, showPlots)
				% Iterate over channels and extract bead data for requested channels
				% in order that they were given, find initial peaks for each channel
				% independently so we can figure out how many bead populations to fit.
				
				extractedBeadData = [];
				figFits = struct();
				counts = struct();
				centers = struct();
				for ch = 1:numel(channels)

					% Extract channel data
					channelData = [];
					for i = 1:length(beadData)
						channelData = [channelData; beadData(i).(channels{ch}).raw]; %#ok<AGROW>
					end
					extractedBeadData(:, ch) = channelData; %#ok<AGROW>
					numBins = round(length(channelData) / 100);

					% Find points greater than 0 since we need to log transform for fitting
					if ~exist('valid', 'var')
						valid = (channelData > 0);
					else
						valid = valid & (channelData > 0);
					end

					% Create histogram w/ biexp data so that peaks are easier to find
					if showPlots
						figFits.(channels{ch}) = figure(); 
						subplot(2, 1, 1)
					end
					[counts.(channels{ch}), centers.(channels{ch})] = Plotting.biexhist( ...
								channelData, numBins, showPlots);
				end
				
				extractedBeadData = extractedBeadData(valid, :);
			end
			
			
			function means = fitGaussians(extractedBeadData)
				
				extrBeadDataBiex = Transforms.lin2logicle(extractedBeadData);
				
				% Have the user draw a set of lines where each point corresponds
				% with the center of a peak
				[~, sortIdx] = sort(max(extractedBeadData, [], 1));
				figPreGauss = figure(); axG = gca(); %#ok<NASGU>
				Plotting.densityplot(axG, ...
						extrBeadDataBiex(:, sortIdx(end)), ...
						extrBeadDataBiex(:, sortIdx(end-1)), ...
						5000, 'hist', ColorMap('parula'));
				title('Draw a line connecting the peaks (in any order)')
				h = impoly();
				position = wait(h);
				numPops = size(position, 1);
				initPops = zeros(length(extractedBeadData), 1);
				for bd = 1:size(extractedBeadData, 1)
					distsX = extrBeadDataBiex(bd, sortIdx(end)) - position(:, 1);
					distsY = extrBeadDataBiex(bd, sortIdx(end-1)) - position(:, 2);
					
					dists = sqrt(distsX.^2 + distsY.^2);
					[~, initPops(bd)] = min(dists);
				end
% 				initPops(1:30)
% 				for ip = 1:numPops, mean(extractedBeadData(initPops == ip, :)), end
				
				% Fit multi-dim gaussian mixture model
				try
					% Fit to biexponential data because data is lognormal, but
					% the low-value data dominates too much if using a typical
					% log transformation. SharedCovariance is true because the
					% covariance is roughly equal in each direction. Could
					% probably also set covariance to diagonal since there is 
					% very little correlation between fluorescence in each peak.
					gmm = fitgmdist(extrBeadDataBiex, numPops, ...
									'Replicates', 1, 'Start', initPops);%, ...
									%'SharedCovariance', true);
				catch ME
					% Ill-conditioned covariance error can come up, which we
					% want to skip, not kill the program
					fprintf(2, '%s\n', ME.message);
					fprintf(2, '--> Error obtained with numPeaks = %d\n', numPops);
					
					means = nan;
					return % Break out
				end
				
				% Sort populations
				means = gmm.mu;
% 				means = (Transforms.logicle2lin(gmm.mu));
				[~, sortIdx] = sort(means(:, 1));
				means = means(sortIdx, :);
			end
			
			
			function [means] = fitFindPeaks(extractedBeadData)
				
				error('Not implemented!\n')
				
			end
						
			
			function filteredMeans = filterMeans(means, counts, centers, thresh)
				% Filters the means based on the max histogram peak 
				
				goodMeans = true(size(means, 1), 1);
				if ~ismember('noFilter', options)
					for ch = 1:numel(channels)
						% Check each channel to see if the mean corresponds with a
						% noise peak, rather than a real peak
						for m = 1:size(means, 1)
							[~, meanLoc] = min(abs(means(m, ch) - centers.(channels{ch})));
							if (counts.(channels{ch})(meanLoc) < (max(counts.(channels{ch})) / thresh))
								goodMeans(m) = false;
							end
						end
					end
				end
				fprintf(1, 'Filtered Means:\n');
				disp(means(goodMeans, :))
				filteredMeans = log10(Transforms.logicle2lin(means(goodMeans, :)));
			end
		end
		
		
		
		function meflFits = calibrateMEFL(controlData, channels, dataType, saveName, showPlots)
			% Computes conversions for MEF units to MEFL units using two-color controls 
			%
			%	Fits are calulated by simple linear regression - only the slope
			%	is recorded, but non-zero intercepts are allowed.
			%
			%	Inputs
			%		controlData		<struct> A standard struct containing only 
			%						two-color control data. 
			%
			%		channels		<cell, char> Cell array of channel names to convert 
			%						to MEF. Must coincide in the same order as the
			%						two-color controls in controlData.
			%						 * Also accepts a char for a single channel
			%
			%		dataType		<char> The dataType to use for calculations
			%						 - This should be the post-MEF conversion 
			%						   compensated dataType (eg 'mComp')
			%
			%		saveName		<char> A name for saving the fits in the
			%						current directory in the 'Calibration' folder 
			%						created by Tranforms.calibrateMEF()
			%
			%		showPlots		(optional) <logical> Flag to show fitting plots
			%
			%	Outputs
			%		meflFits		<struct> A struct array where each field is
			%						the channel name of a converted channel and the 
			%						value is the conversion factor to get MEFL units. 
			
			checkInputs_calibrateMEFL();
			
			% Check if conversions already exist
			beadDir = ['Calibration', filesep];
			if ~exist(beadDir, 'file')
				error('Bead directory not found! It should be set up during MEF calibration')
			end
			meflFname = [beadDir, saveName, '_MEFL_Conversions.mat'];
			if exist(meflFname, 'file')
				fprintf(1, 'Loading pre-computed MEFL conversions\n');
				load(meflFname)
				return
			end
			
			% Compute conversions
			meflFits = struct();
			for ch = 1:numel(channels)
				currChan = channels{ch};
				if strcmpi(currChan, 'FITC_A') 
					% Skip FITC channel since it is the reference channel
					meflFits.(currChan) = 1;
				else
					currChanMEF = Transforms.CHANNEL_MAP.(currChan);
					
					xdata = controlData(ch).(currChan).(dataType);
					ydata = controlData(ch).FITC_A.(dataType);
					pos = ((xdata > 0) & (ydata > 0));
					xdataLog = log10(xdata(pos));
					ydataLog = log10(ydata(pos));
					ss_tot = sum((ydataLog - mean(ydataLog)).^2);
					
					% Fit in log-space so that highly transfected cells don't
					% completely dominate
					opt = optimoptions('lsqcurvefit', 'Display', 'off');
					[fit, ss_res] = lsqcurvefit(@(p, x) p(2) + x, [1, 1], ...
							xdataLog, ydataLog, [], [], opt);
					slope = 10.^fit(2);
					rsq = 1 - (ss_res / ss_tot);
					
					% No longer using linear fitting due to overstated influence
					% of highly transfected cells
% 					[rval, slope, ~] = regression(xdata', ydata'); % Convert column vectors to rows
% 					rsq = rval.^2;

					meflFits.(currChan) = slope;
					
					if (exist('showPlots', 'var') && showPlots)
						figMEFL = figure(); 
						ax = gca(); hold(ax, 'on');
						xrange = logspace(-2, 7, 50);
						
% 						plot(ax, Transforms.lin2logicleMEF(xdata), ...
% 							 Transforms.lin2logicleMEF(ydata), ...
% 							 '.', 'MarkerSize', 2)
% 						plot(ax, Transforms.lin2logicleMEF(xrange), ...
% 							 Transforms.lin2logicleMEF(xrange * slope), ...
% 							 '-', 'LineWidth', 4)
						 
						plot(ax, xdata, ydata, ...
							 '.', 'MarkerSize', 2)
						plot(ax, xrange, xrange * slope, ...
							 '-', 'LineWidth', 4)
						
% 						Plotting.biexpAxesMEF(ax);
						ax.YScale = 'log';
						ax.XScale = 'log';
						title(sprintf('Scale Factor: %.2f | R^2: %.3f', ...
								slope, rsq), 'fontsize', 14)
						ylabel('MEFL')
						xlabel(currChanMEF);
						
						figFname = [beadDir, saveName, '_', currChanMEF, '_MEFL_Conversion.fig'];
						saveas(figMEFL, figFname)
					end
				end
			end
			save(meflFname, 'meflFits')
			
			
			% --- Helper Functions --- %
			
			
			function checkInputs_calibrateMEFL()
				validateattributes(controlData, {'struct'}, {}, mfilename, 'controlData', 1);
				validateattributes(channels, {'cell', 'char'}, {}, mfilename, 'channels', 2);
				if ischar(channels), channels = {channels}; end % Convert to cell for simplicity
				assert(numel(channels) == numel(controlData), 'Number of channels and controls must be equal!\n');
				
				validChannels = fieldnames(controlData(1));
				badChannels = setdiff(channels, validChannels);
				assert(isempty(badChannels), 'Channel not in controlData: %s\n', badChannels{:});
				
				validatestring(dataType, fieldnames(controlData(1).(channels{1})), mfilename, 'dataType', 3);
				validateattributes(saveName, {'char'}, {}, mfilename, 'saveName', 4);
			end
		end
		
		
		
        function data = fcs2ABC(data, channelFits, dataType)
			% Converts fcs data to ABC units
            %
			%	data = fcs2ABC(data, channelFits, dataType)
			%
            %   Inputs
            %
			%       data            A standard data struct. The channels must match 
			%						those in channelFits. 
            %
            %       channelFits     (Output from calibrateABC())
            %                       An Nx2 table of slopes/intercepts for linear conversion of 
            %                       N channels to ABC equivalents. Rows are labeled with 
            %                       channel names and columns are labeled with slope / zero
            %
            %       dataType        The population to convert to ABCs (eg 'scComp' or 'tf')
            %                        - must be in data!
			%
			%
			%	Outputs
			%
			%		data			Data struct with ABC units added to all converted 
			%						channels in a new field called 'abc'.
			%
			% Written by
			% Ross Jones
			% jonesr18@mit.edu
			% 2016-04-23
            
            % Process inputs
			channels = checkInputs_fcs2ABC(data, channelFits, dataType);
                        
			for i = 1:numel(data) %#ok<ALIGN>
                
                for ch = 1:numel(channels)
                
                    channel = channels{ch};
                    chData = data(i).(channel).(dataType);
                    
                    % Remove negative values -- might be better to allow isolation post-proc
                    if any(chData < 0)
                        chData = chData(chData >= 0);
                    end
                    
                    % Do conversion
                    m = channelFits(channel, :).slope;
                    b = channelFits(channel, :).zero;
                    data(i).(channel).abc = 10^b .* chData.^m;                    
                end
            end
			
			
			% --- Helper Functions --- %
			
            
			function channels = checkInputs_fcs2ABC(data, channelFits, dataType)
				
				validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                validateattributes(channelFits, {'table'}, {}, mfilename, 'channelFits', 2);
                
                tableChannels = channelFits.Properties.RowNames;
                dataChannels = data(1).chanNames;
                if ~isempty(setdiff(tableChannels, dataChannels))
                    error('Channels in channelFits do not match channels in data!')
                end
                
                channels = tableChannels;
                validatestring(dataType, fieldnames(data(1).(channels{1})), mfilename, 'dataType', 3);
			end

        end
        
        
        function channelFits = calibrateABCfromMedians(beadsFilenames, channel, hasBlank)
            % calibrateABCfromMedians creates antibody binding capacity (ABC) fits for an 
            % .fcs file of a particular bead sample supplied by the user using medians of 
            % each file. Note this is for Quantum Simply Cellular (R) beads, not rainbow 
            % beads (see calibrateMEF()). 
            %
            % Calibration plots and the table with ABC fits are saved in pwd().
            %
            %   Inputs: 
            %   
            %       beadFilename        A cell array of .fcs filenames which contain Ab-stained 
            %                           bead data. Can be given in any order (median is sorted).
            %
            %       channel             The channel name to calibrate bead data
            %
            %       hasBlank            (Optional) Boolean to indicate if the samples include the
            %                           blank population (ignore lowest peak if so). Default = false
            %
            %   Outputs:
            %
            %       channelFit          An Nx2 table of slopes/intercepts for linear conversion of 
            %                           N channels to ABC equivalents. Rows are labeled with 
            %                           channel names and columns are labeled with slope / zero
            %       
            %       medians             A 1xN array of medians values in order the files are given
            %
            % Written by Ross Jones
            % 2016-06-01
            % MIT Weiss Lab
            % jonesr18@mit.edu
            
            % Import data
            beadData = FlowAnalysis.openFiles(beadsFilenames{:});
            
            % Initialize as an array, convert to a table later
            channelFits = zeros(1, 2);
            
            % Calculate medians
            beadMedians = zeros(1, numel(beadData));
            for i = 1:numel(beadData)
                % Check if unimodal (if not, there are many unstained beads)
%                 sortedBeads = sort(beadData(i).(channel).raw);
%                 logSortedBeads = Transforms.lin2logicle(sortedBeads);
%                 [dip, p_value] = hartigansDipSignifTest(logSortedBeads, 5000);
%                 if (p_value < 0.05)
%                     warning('Significant non-unimodal behavior detected: Dip = %.2f, P = %2f', ...
%                         dip, p_value);
%                     
%                     % If multi-modal, compute mean of higher population (stained beads) with
%                     % gaussian mixture model
%                     figure()
%                     Plotting.biexhist(sortedBeads);
%                     gmm = fitgmdist(logSortedBeads, 2, 'Start', 'plus', 'Replicates', 5);
%                     beadMedians(i) = Transforms.logicle2lin(max(gmm.mu));
%                 else
                    beadMedians(i) = median(beadData(i).(channel).raw);
%                 end
            end
            beadMediansSorted = sort(Transforms.lin2logicle(beadMedians));
            
            % Extract channel data and plot as histogram
            figPeaks = figure();
            ax1 = subplot(1, 2, 1);
            hold(ax1, 'on');
            for i = 1:length(beadData)
                channelData = beadData(i).(channel).raw; 
                numBeads = length(channelData);
                numBins = numBeads / 100;
				% Create histogram w/ biexp data so that peaks are easier to find
                [binCounts, binCenters] = Plotting.biexhist(channelData, numBins);
                
                % Match closest bin to medians
                [~, binIdx] = min(abs(binCenters - Transforms.lin2logicle(beadMedians(i))));
                peak = binCounts(binIdx);
                plot(ax1, Transforms.lin2logicle(beadMedians(i)), peak, '*r')
            end         
            title('Bead Peak Medians')
            ylabel('COUNT')
            xlabel(strrep(channel, '_', '-'))
            ax1.FontSize = 14;
            hold(ax1, 'off')
            
            % Assume we have the P highest peaks and an unstained population
            ABC = [5899, 68458, 281582, 621651];
            if (exist('hasBlank', 'var') && hasBlank)
                P = min(length(beadMediansSorted) - 1, numel(ABC));
            else
                P = min(length(beadMediansSorted), numel(ABC));
            end
            ABC = ABC(end - P + 1 : end);
            peakIntensity = Transforms.logicle2lin(beadMediansSorted); % Revert to linear for fitting
            peakIntensity = peakIntensity(end - P + 1 : end);
            
            % Calculate fit
            fit = polyfit(log10(peakIntensity), log10(ABC), 1);
            
            ax2 = subplot(1, 2, 2);
            plot(ax2, ...
                peakIntensity, ABC, 'ob', ...
                peakIntensity, 10.^polyval(fit, log10(peakIntensity)), 'r-')
            ax2.XScale = 'log';
            ax2.YScale = 'log';
            ax2.XLim = [1e-1, 1e5];
            ax2.YLim = [1e3, 1e7];
            ax2.FontSize = 14;
            
            % Calculate R^2 value
            yresid = ABC - 10.^polyval(fit, log10(peakIntensity));
            SSresid = sum(yresid.^2);
            SStotal = (length(ABC) - 1) * var(ABC);
            rsq = 1 - SSresid / SStotal;
            title(sprintf('R^2 = %.3f', rsq))
            xlabel('Fluorescence')
            ylabel('ABC')
            hold(ax2, 'off')

            % Assign fit data
            channelFits(1, :) = fit;

            fprintf(1, 'Finished fitting %s. # Peaks = %d\n', channel, length(peaks));
            savefig(figPeaks, [channel '_ABC_Fit.fig'])

            % Convert fits to a table
            channelFits = array2table(channelFits);
            channelFits.Properties.VariableNames = {'slope', 'zero'};
            channelFits.Properties.RowNames = {channel};

            % Save channelFits
            save([channel, '_ABC_Fit.mat'], 'channelFits')
        end
        
        
        function [channelFits, locs, peaks] = calibrateABC(beadsFilename, channel, largePeaks, hasBlank, locs, peaks)
            % calibrateABC creates antibody binding capacity (ABC) fits for an .fcs file 
            % of a particular bead sample supplid by the user. Note this is for Quantum 
            % Simply Cellular (R) beads, not rainbow beads (see calibrateMEF()). 
            %
            % Calibration plots and the table with ABC fits are saved in pwd().
            %
            %   Inputs: 
            %   
            %       beadFilename        An .fcs filename which contains Ab-stained bead data
            %                           If a cell array of strings, the data from each file 
            %                           will be concatenated together
            %
            %       channel             The channel name to calibrate bead data
            %
            %       largePeaks          (Optional) Boolean to indicate whether to rule out peaks
            %                           smaller than 1/10 max(peaks) on the first peak finding 
            %                           run. Default = true
            %
            %       hasBlank            (Optional) Boolean to indicate if the samples include the
            %                           blank population (ignore lowest peak if so). Default = false
            %
            %       locs                (Optional) 1xN array of peak locations in ascending order
            %
            %       peaks                 (Optional) 1xN array of peak heights
            %
            %   Outputs:
            %
            %       channelFit          An Nx2 table of slopes/intercepts for linear conversion of 
            %                           N channels to ABC equivalents. Rows are labeled with 
            %                           channel names and columns are labeled with slope / zero
            %
            %       locs                A 1xN array of peak locations in ascending order
            %
            %       peaks                 A 1xN array of peak heights
            %
            % Written by Ross Jones
            % 2016-04-23
            % MIT Weiss Lab
            % jonesr18@mit.edu
            % 
            % Updated 2016-06-01
            
            % Iterate over all given filenames and concatenate data together. 
%             assert(iscell(beadsFilenames), 'beadsFilename must be a cell array of strings!');
            if ~iscell(beadsFilename)
                beadsFilename = {beadsFilename};
            end
            beadData = FlowAnalysis.openFiles(beadsFilename{:});
            if ~exist('largePeaks', 'var')
                largePeaks = true;
            end
            
            % Initialize as an array, convert to a table later
            channelFits = zeros(1, 2);

            % Extract channel data
            channelData = [];
            for i = 1:length(beadData)
                channelData = [channelData; beadData(i).(channel).raw]; %#ok<AGROW>
            end
            numBeads = length(channelData);
            
            % Plot figure to show fitting
            figPeaks = figure('Position', [500 500 700 350]);
            subplot(1, 2, 1)
            numBins = round(numBeads / 100);
            [counts, centers] = Plotting.biexhist(channelData, numBins);
            hold on

            % Find highest peaks
            if ~(exist('locs', 'var') && exist('peaks', 'var'))
                [peaks, locs] = findpeaks(counts);
                if (largePeaks)
                    maxpeak = max(peaks);
                    [peaks, locs] = findpeaks(counts, 'MinPeakHeight', maxpeak / 10);
                end
            end
            
            % Plot peaks on data
            subplot(1, 2, 1)
            hold on
            plot(centers(locs), peaks, '*r')
            title(strrep(channel, '_', '-'))
            ylabel('COUNT')
            
            % Assume we have the P highest peaks and an unstained population
            ABC = [68458, 281582, 621651];
%             ABC = [5899, 68458, 281582, 621651];

            if (exist('hasBlank', 'var') && hasBlank)
                P = min(length(locs) - 1, numel(ABC));
            else
                P = min(length(locs), numel(ABC));
            end
            ABC = ABC(end-P+1:end);
            peakIntensity = Transforms.logicle2lin(centers(locs));
            peakIntensity = peakIntensity(end-P+1:end);
            
            % Calculate fit
            r = 2^18;   % resolution
            n = 5;      % log decades
            peakRelChannel = r / n * log10(peakIntensity);
            fit = polyfit(peakRelChannel, log10(ABC), 1);
            
            subplot(1,2,2)
            semilogy(peakRelChannel, ABC, 'ob')
            hold on
            semilogy(peakRelChannel, 10.^polyval(fit, peakRelChannel), 'r-')
            
            % Calculate R^2 value
            yresid = ABC - 10.^polyval(fit, peakRelChannel);
            SSresid = sum(yresid.^2);
            SStotal = (length(ABC) - 1) * var(ABC);
            rsq = 1 - SSresid / SStotal;
            title(sprintf('R^2 = %.3f', rsq))
            axis([0 2^18 1 1e6])
            xlabel('REL CH #')
            ylabel('ABC')
            hold off

            % Assign fit data
            channelFits(1, :) = fit;

            fprintf(1, 'Finished fitting %s. # Peaks = %d\n', channel, length(peaks));
            savefig(figPeaks, [channel '_ABC_Fit.fig'])

            % Convert fits to a table
            channelFits = array2table(channelFits);
            channelFits.Properties.VariableNames = {'slope', 'zero'};
            channelFits.Properties.RowNames = {channel};

            % Save channelFits
            save([channel, '_ABC_Fit.mat'], 'channelFits')
        end
        
    end
    
end