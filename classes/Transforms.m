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
    % Written/Compiled by
	% Ross Jones
	% jonesr18@mit.edu
    % Weiss Lab, MIT
	
	properties (Constant)
		
		MEF_CONVERSION_FACTOR = 1e3;
		
		% Unique names - has error when you select more peaks than the one w/
		% the least number of peaks!
% 		CHANNEL_MAP = struct( ...
% 			'BUV_396_A',		'MEPB', ...		% Not exactly right but the closest one
% 			'Cascade_Blue_A',	'MEPB', ...
%			'Pacific_Blue_A',	'MEPB', ...
%			'AmCyan_A',			'MEAmC', ...
% 			'FITC_A',			'MEFL', ...
% 			'PE_A',				'MEPE', ...
%			'PE_YG_A'			'MEPE', ...		% Alt name for Koch LSRII-HTS2
% 			'PE_Texas_Red_A',	'MEPTR', ...
% 			'PE_TxRed_YG_A',	'MEPTR', ...	% Alt name for Koch LSRII-HTS2
% 			'PE_Cy5_5_A',		'MECY', ...
% 			'PE_Cy7_A',			'MEPCY7', ...
% 			'APC_A',			'MEAP', ...
% 			'APC_Cy7_A',		'MEAPCY7');
		
		% All MEFLs - ignore peak # error, calculation should come out the same
		CHANNEL_MAP = struct( ...
			'BUV_396_A',		'MEFL', ...
			'Cascade_Blue_A',	'MEFL', ...
			'Pacific_Blue_A',	'MEFL', ...
			'AmCyan_A',			'MEFL', ...
			'FITC_A',			'MEFL', ...
			'PE_A',				'MEFL', ...
			'PE_YG_A',			'MEFL', ...		% Alt name for Koch LSRII-HTS2
			'PE_Texas_Red_A',	'MEFL', ...
			'PE_TxRed_YG_A',	'MEFL', ...		% Alt name for Koch LSRII-HTS2
			'PE_Cy5_5_A',		'MEFL', ...
			'PE_Cy7_A',			'MEFL', ...
			'APC_A',			'MEFL', ...
			'APC_Cy7_A',		'MEFL');
	end
	
	methods (Static)
		
		function [biexpData, params] = lin2logicle(S, doMEF, params)
			% Converts input linear data (S) to a logicle scale 
			% 
			%	biexpData = Transforms.lin2logicle(linspace(0, 1000));
			%
			%	Inputs
			%		S			<numeric> Linear-scale input values
			%
			%		doMEF		<logical> (optional) If true, the method will scale
			%					the transformation for MEF (bead unit) conversions,
			%					rather than raw machine-based fluorescent units
			%
			%		params		<struct> (optional) Allows the user to specify
			%					the exact logicle transform parameters to use
			%						Accepted input fields (see below):
			%							T, M, r, MEF
			%
			%	Outputs
			%		biexpData	<numeric> Logicle-transformed data
			%
			%		params		<struct> The parameters used for conversion
			%
			%	Notes
			%		The logicle conversion function is described by Parks, et al. 
			%		"A New Logicle Display Method Avoids Deceptive Effects of 
			%		Logarithmic Scaling for Low Signals and Compensated Data" 
			%   
			%		This method implements the inverse of Equation (5):
			%	
			%			S = T * 10^(M-W) * (10^(X-W) - p^2 * 10^-((X-W)/p) + p^2 - 1)
			%				(for X >= W)
			%		
			%		Where (see Figure 2):
			%			S := linear 'raw' values
			%			X := logicle-scale values
			%			T := "Top" of scale data value (ie the highest data
			%				 value expected from the machine, typically 2^18)
			%			M := Total plot width in asmymptototic decades
			%				 (defualt = 4.5)
			%			W := Width of linearization, computed as such:
			%				 W = (M - log10(T/abs(r))) / 2		| Eqn (6)
			%				 This determines the value of lin2logicle(0)
			%			r := Reference point for negative data range
			%				 (default = -150)
			%
			%		Since Equation 5 cannot be explicitly solved for X, a 
			%		spline is fitted to generic X data to estimate the logicle
			%		values of S
			%
			%		For MEF (bead unit) data, the default is to scale the data
			%		values by Transforms.MEF_CONVERSION_FACTOR prior to doing
			%		the logicle transformation. This parameter can be
			%		overwritten with the optional params.MEF input.
			%
			%
			% Written By
			% Breanna DiAndreth
			% bstillo@mit.edu
			% Weiss Lab, MIT
			%
			% Update Log:
			%	2018-02-10 (Ross):	Merged this w/ lin2logicleMEF to do MEF 
			%						conversion optionally, also added optional
			%						logicle parameters input
			
			% Check inputs
			zCheckInputs_lin2logicle();
			
			% Converting logicle to linear is solveable, but the reverse is not, 
			% so we fit a spline to generic logicle-lin data and then project
			% the given data onto that spline
			X = linspace(0, params.M, 1000);				% Generic X values
			SX = Transforms.logicle2lin(X, doMEF, params);	% S values from X
			p = spline(SX / params.MEF, X);					% Spline fit
			
			% Apply spline fit to the actual given S values
			biexpData = ppval(p, S / params.MEF);						
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_lin2logicle()
				
				% Check input arguments
				validateattributes(S, {'numeric'}, {}, mfilename, 'S', 1);
				
				doMEF = (exist('doMEF', 'var') && all(logical(doMEF)));
				
				if exist('params', 'var')
					validateattributes(params, {'struct'}, {}, mfilename, 'params', 3);
				else
					params = struct();
				end
				params = Transforms.checkLogicleParams(doMEF, params);
			end
		end
		
		
		function [logicleData, params] = lin2logicle2(linearData, scale, params)
			% Converts input linear data to a logicle scale 
			% 
			%	biexpData = Transforms.lin2logicle2(linearData, scale, params);
			%
			%	Inputs
			%		linearData	<numeric> Linear-scale input values
			%
			%		scale		<numeric> A scaling factor for 'stretching' the
			%					transform to cover scaled values (e.g. MEFLs)
			%
			%		params		<struct> (optional) Allows the user to specify
			%					the exact logicle transform parameters to use
			%						Accepted input fields (see below):
			%							T, M, r
			%
			%	Outputs
			%		logicleData	<numeric> Logicle-transformed data
			%
			%		params		<struct> The parameters used for conversion
			%
			%	Notes
			%		The logicle conversion function is described by Parks, et al. 
			%		"A New Logicle Display Method Avoids Deceptive Effects of 
			%		Logarithmic Scaling for Low Signals and Compensated Data" 
			%   
			%		This method implements the inverse of Equation (5):
			%	
			%			S = T * 10^(M-W) * (10^(X-W) - p^2 * 10^-((X-W)/p) + p^2 - 1)
			%				(for X >= W)
			%		
			%		Where (see Figure 2):
			%			S := linear 'raw' values
			%			X := logicle-scale values
			%			T := "Top" of scale data value (ie the highest data
			%				 value expected from the machine, typically 2^18)
			%			M := Total plot width in asmymptototic decades
			%				 (defualt = 4.5)
			%			W := Width of linearization, computed as such:
			%				 W = (M - log10(T/abs(r))) / 2		| Eqn (6)
			%				 This determines the value of lin2logicle(0)
			%			r := Reference point for negative data range
			%				 (default = -150)
			%
			%		Since Equation 5 cannot be explicitly solved for X, a 
			%		spline is fitted to generic X data to estimate the logicle
			%		values of S
			%
			%
			% Written By
			% Breanna DiAndreth
			% bstillo@mit.edu
			% Weiss Lab, MIT
			%
			% Update Log:
			%	2018-02-10 (Ross):	Merged this w/ lin2logicleMEF to do MEF 
			%						conversion optionally, also added optional
			%						logicle parameters input
			%	2020-02-26 (Ross):	Updated 'scale' input to match biexpAxes2
			
			% Check inputs
			zCheckInputs_lin2logicle2();
			
			% Converting logicle to linear is solveable, but the reverse is not, 
			% so we fit a spline to generic logicle-lin data and then project
			% the given data onto that spline
			X = linspace(0, params.M, 1000);				% Generic X values
			SX = Transforms.logicle2lin2(X, scale, params);	% S values from X
			p = spline(SX ./ scale, X);					% Spline fit
			
			% Apply spline fit to the actual given S values
			logicleData = ppval(p, linearData ./ scale);						
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_lin2logicle2()
				
				% Check input arguments
				validateattributes(linearData, {'numeric'}, {}, mfilename, 'linearData', 1);
				
				if exist('scale', 'var')
					validateattributes(scale, {'numeric'}, {'scalar'}, mfilename, 'scale', 2);
				else
					scale = 1;
				end
				
				if exist('params', 'var')
					validateattributes(params, {'struct'}, {}, mfilename, 'params', 3);
				else
					params = struct();
				end
				
				params = Transforms.checkLogicleParams2(params);
			end
		end
		
		
		function [linData, params] = logicle2lin(X, doMEF, params)
			% Converts input logicle data (X) to a linear scale 
			% 
			%	linData = Transforms.logicle2lin(linspace(0, 4.5));
			%
			%	Inputs
			%		X			<numeric> Logicle-scale input values
			%
			%		doMEF		<logical> (optional) If true, the method will scale
			%					the transformation for MEF (bead unit) conversions,
			%					rather than raw machine-based fluorescent units
			%
			%		params		<struct> (optional) Allows the user to specify
			%					the exact logicle transform parameters to use
			%						Accepted input fields (see below):
			%							T, M, r, MEF
			%
			%	Outputs
			%		linData		<numeric> Linear-transformed data
			%
			%		params		<struct> The parameters used for conversion
			%
			%	Notes
			%		The logicle conversion function is described by Parks, et al. 
			%		"A New Logicle Display Method Avoids Deceptive Effects of 
			%		Logarithmic Scaling for Low Signals and Compensated Data" 
			%   
			%		This method implements Equation (5):
			%	
			%			S = T * 10^(M-W) * (10^(X-W) - p^2 * 10^-((X-W)/p) + p^2 - 1)
			%				(for X >= W)
			%		
			%		Where (see Figure 2):
			%			S := linear 'raw' values
			%			X := logicle-scale values
			%			T := "Top" of scale data value (ie the highest data
			%				 value expected from the machine, typically 2^18)
			%			M := Total plot width in asmymptototic decades
			%				 (defualt = 4.5)
			%			W := Width of linearization, computed as such:
			%				 W = (M - log10(T/abs(r))) / 2		| Eqn (6)
			%				 This determines the value of lin2logicle(0)
			%			r := Reference point for negative data range
			%				 (default = -150)
			%
			%		For MEF (bead unit) data, the default is to scale the data
			%		values by Transforms.MEF_CONVERSION_FACTOR after returning
			%		from the logicle transformation. This parameter can be
			%		overwritten with the optional params.MEF input.
			%
			% Written By
			% Breanna DiAndreth
			% bstillo@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%   2015-02-09 (Ross):	Added comments and intuitive variable names
			%						Major speed up by using algebraic conversion 
            %						instead of symbolic toolbox
            %   2016-03-28 (Ross):	Made display breadth smaller w/ Bre's adjustment 
			%						of r, which increases the relative size of the
            %						area between 10^-2 and 10^2
			
			zCheckInputs_logicle2lin()
			
			% This gives a range for linearization around zero (W).
			W = (params.M - log10(params.T / abs(params.r))) / 2;

			% p and W are considered one parameter, p is introduced by the authors
			% for compactness. Here we find p that solves their equivalence:
			%   w = W * ln(10) = 2 * p * ln(p) / (p + 1)
			% Solved by WolframAlpha:
			w = W * log(10);
			p = w / (2 * lambertw(w / 2 * exp(-w / 2)));
			
			% Find where the given data is below zero (thus invalid).
			logX = (X >= W);
			
			% Compute and return the final linearized vector 
			linData = params.MEF .* (toLin(X, params) .* logX ...
						  - toLin(2 * W - X, params) .* (1 - logX));
			
			
			% --- Helper Functions --- %
			
			
			function S = toLin(X, params)
				% Does the final conversion		| Eqn (5)
				
				S = params.T .* 10^-(params.M - W) .* ...
					(10.^(X - W) - p^2 .* 10.^(-(X - W) ./ p) + p^2 - 1);
			end
						
			
			function zCheckInputs_logicle2lin()
				
				% Check input arguments
				validateattributes(X, {'numeric'}, {}, mfilename, 'X', 1);
				
				doMEF = (exist('doMEF', 'var') && all(logical(doMEF)));
				
				if exist('params', 'var')
					validateattributes(params, {'struct'}, {}, mfilename, 'params', 3);
				else
					params = struct();
				end
				params = Transforms.checkLogicleParams(doMEF, params);

			end
		end
		
		
		function [linearData, params] = logicle2lin2(logicleData, scale, params)
			% Converts input logicle data to a linear scale 
			% 
			%	linearData = Transforms.logicle2lin2(logicleData, scale, params);
			%
			%	Inputs
			%		logicleData <numeric> Logicle-scale input values
			%
			%		scale		<numeric> A scaling factor for 'stretching' the
			%					transform to cover scaled values (e.g. MEFLs)
			%
			%		params		<struct> (optional) Allows the user to specify
			%					the exact logicle transform parameters to use
			%						Accepted input fields (see below):
			%							T, M, r
			%
			%	Outputs
			%		linearData	<numeric> Linear-transformed data
			%
			%		params		<struct> The parameters used for conversion
			%
			%	Notes
			%		The logicle conversion function is described by Parks, et al. 
			%		"A New Logicle Display Method Avoids Deceptive Effects of 
			%		Logarithmic Scaling for Low Signals and Compensated Data" 
			%   
			%		This method implements Equation (5):
			%	
			%			S = T * 10^(M-W) * (10^(X-W) - p^2 * 10^-((X-W)/p) + p^2 - 1)
			%				(for X >= W)
			%		
			%		Where (see Figure 2):
			%			S := linear 'raw' values
			%			X := logicle-scale values
			%			T := "Top" of scale data value (ie the highest data
			%				 value expected from the machine, typically 2^18)
			%			M := Total plot width in asmymptototic decades
			%				 (defualt = 4.5)
			%			W := Width of linearization, computed as such:
			%				 W = (M - log10(T/abs(r))) / 2		| Eqn (6)
			%				 This determines the value of lin2logicle(0)
			%			r := Reference point for negative data range
			%				 (default = -150)
			%
			%		For MEF (bead unit) data, the default is to scale the data
			%		values by Transforms.MEF_CONVERSION_FACTOR after returning
			%		from the logicle transformation. This parameter can be
			%		overwritten with the optional params.MEF input.
			%
			% Written By
			% Breanna DiAndreth
			% bstillo@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%   2015-02-09 (Ross):	Added comments and intuitive variable names
			%						Major speed up by using algebraic conversion 
            %						instead of symbolic toolbox
            %   2016-03-28 (Ross):	Made display breadth smaller w/ Bre's adjustment 
			%						of r, which increases the relative size of the
            %						area between 10^-2 and 10^2
			%	2020-02-26 (Ross):	Updated 'scale' input to match biexpAxes2
			
			zCheckInputs_logicle2lin2()
			
			% This gives a range for linearization around zero (W).
			W = (params.M - log10(params.T / abs(params.r))) / 2;
			
			% p and W are considered one parameter, p is introduced by the authors
			% for compactness. Here we find p that solves their equivalence:
			%   w = W * ln(10) = 2 * p * ln(p) / (p + 1)
			% Solved by WolframAlpha:
			w = W * log(10);
			p = w / (2 * lambertw(w / 2 * exp(-w / 2)));
			
			% Find where the given data is above (thus valid).
			posData = (logicleData >= W);
			
			% Compute and return the final linearized vector 
			linearData = scale .* (toLin(logicleData, params) .* posData ...
						  - toLin(2 * W - logicleData, params) .* (1 - posData));
			
			
			% --- Helper Functions --- %
			
			
			function S = toLin(X, params)
				% Does the final conversion		| Eqn (5)
				
				S = params.T .* 10^-(params.M - W) .* ...
					(10.^(X - W) - p^2 .* 10.^(-(X - W) ./ p) + p^2 - 1);
			end
						
			
			function zCheckInputs_logicle2lin2()
				
				% Check input arguments
				validateattributes(logicleData, {'numeric'}, {}, mfilename, 'logicleData', 1);
				
				if exist('scale', 'var')
					validateattributes(scale, {'numeric'}, {'scalar'}, mfilename, 'scale', 2);
				else
					scale = 1;
				end
				
				if exist('params', 'var')
					validateattributes(params, {'struct'}, {}, mfilename, 'params', 3);
				else
					params = struct();
				end
				
				params = Transforms.checkLogicleParams2(params);

			end
		end
		
		
		function requestedUnits = getBeadUnits(channels)
			% Returns the MEF-equivalent unit names for the given channels
			%
			%	requestedUnits = getBeadUnits(channels)
			%
			%	Inputs
			%		channels		<cell, char> The channels to extract bead
			%						values for (eg 'Pacific_Blue_A')
			%
			%	Outputs
			%		requestedUnits	<cell> A cell array of string MEF unit names
			%						in the same order as 'channels'
			% 
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			% 
			
			% Check inputs
			validateattributes(channels, {'char', 'cell'}, {}, mfilename, 'channels', 1);
			if ischar(channels), channels = {channels}; end % For simplicity
			badChannels = setdiff(channels, fieldnames(Transforms.CHANNEL_MAP));
			assert(isempty(badChannels), ...
				'Channel not valid: %s', badChannels{:});
			
			% Convert channel names to MEF units
			requestedUnits = cell(1, numel(channels));
			for i = 1:numel(channels)
				requestedUnits{i} = Transforms.CHANNEL_MAP.(channels{i}); 
			end
			
			if numel(unique(requestedUnits)) < numel(channels)
				warning('Non-unique bead units detected: check channel names if this is not intentional')
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
			% 
			%	The valid bead lots and the respective fluorophore counts are in
			%	the folder /Spherotech Bead Data/. 
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			% 
						
			% Find files w/ bead data: 'BEADS_<beadType>_<beadLot1_beadLot2...>.txt'
			thisFname = mfilename('fullpath');
			fnameSplit = split(thisFname, filesep);
			folderName = join(fnameSplit(1:end-2), filesep); % Comes out as a 1x1 Cell
			beadValueFiles = dir([folderName{:}, filesep, '**', filesep, 'BEADS*.txt']);
			beadValueFnames = {beadValueFiles.name};
			
			% Extract bead types and lots
			beadValueFnamesSplit = cellfun(@(x) split(x, {'_', '.'}), ...
					beadValueFnames, 'uniformoutput', false);
			beadTypes = cellfun(@(x) x{2}, beadValueFnamesSplit, 'uniformoutput', false);
			validBeadTypes = unique(beadTypes);
			beadLots = cellfun(@(x) x(3:end-1)', beadValueFnamesSplit, 'uniformoutput', false);
			validBeadLots = unique(cat(2, beadLots{:}));
			
			% Check input types
			validatestring(beadType, validBeadTypes, mfilename, 'beadType', 1);
			validatestring(beadLot, validBeadLots, mfilename, 'beadLot', 2);
			validateattributes(channels, {'char', 'cell'}, {}, mfilename, 'channels', 3);
			
			% Check channels are valid
			if ischar(channels), channels = {channels}; end % For simplicity
			badChannels = setdiff(channels, fieldnames(Transforms.CHANNEL_MAP));
			assert(isempty(badChannels), ...
				'Channel not valid: %s', badChannels{:});
			requestedUnits = Transforms.getBeadUnits(channels);
			
			% Find which bead value file to load data from
			whichFile = beadValueFnames{ ...
					contains(beadValueFnames, beadType) & ...
					contains(beadValueFnames, beadLot)};
			beadValuesTable = readtable([beadValueFiles(1).folder, filesep, whichFile], ...
					'delimiter', '\t');
			beadValues = table2struct(beadValuesTable, 'ToScalar', true);
			
			% Collect values into struct
			beadVals = struct();
			for u = requestedUnits
				nonnan = ~isnan(beadValues.(u{:}));
				beadVals.(u{:}) = beadValues.(u{:})(nonnan);
			end
		end

		
		function data = fcs2MEF(data, channelFits, dataType)
			% Converts fcs data to MEF units
			%	For conversion information, bead types/lots, valid channels, and
			%	MEF unit names, see Transforms.getBeadVals().
            %
            %   Inputs
            %
			%       data            A standard data struct. The channels must be present in channelFits. 
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
			% Written By
			% Breanna Stillo
			% bstillo@mit.edu
			% Weiss Lab, MIT
			%
            % Update Log: 
			%	2016-04-23 by Ross Jones
            %       Uses standard struct interface and uses channelFits directly rather than
            %       loading the saveFile
			
            % Process inputs
			channels = zCheckInputs_fcs2MEF();
                        
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
			
            
			function channels = zCheckInputs_fcs2MEF()
				
				validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                validateattributes(channelFits, {'table'}, {}, mfilename, 'channelFits', 2);
                
                channels = setdiff(fieldnames(data(1)), {'gates', 'nObs', 'Time'});
                badChannels = setdiff(channels, channelFits.Properties.VariableNames);
                assert(isempty(badChannels), ...
                    'Channel not in fits: %s\n', badChannels{:});
                
                validatestring(dataType, fieldnames(data(1).(channels{1})), mfilename, 'dataType', 3);
			end

        end
		
        
        function [mefFits, figFits] = calibrateMEF(beads, channels, options)
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
            %       mefFits			<table> An Nx2 table of slopes/intercepts for  
            %                       linear conversion of N channels to MEF equivalents. 
            %                       Rows are labeled with channel names and columns are 
			%						labeled with slope / zero.
            %       figFits         <struct> A struct of handles to all figures
            %						generated by the method. If 'showPlots' is not
			%						included as an optional input, this will return
			%						as empty.
            %
            % Written By
            % Breanna DiAndreth & Ross Jones
            % bstillo@mit.edu
			% jonesr18@mit.edu
            % 
			% Update Log:
			%	2017-09-16 by Ross Jones
			%		Overhauled interface and processing
			
            % Check inputs, initialize data struct
            [beadData, beadVals] = zCheckInputs_calibrateMEF();
			MEF_units = fieldnames(beadVals)'; % Start in same order as channels
			if (numel(MEF_units) == 1 && numel(channels) > 1)
				warning('Adjusting %s units to match all channels', MEF_units{1});
				MEF_units = repmat(MEF_units, 1, numel(channels));
			elseif (numel(MEF_units) < numel(channels))
				error('Number of bead units does not match number of channels!')
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
				[means, figFits.manualPeaks] = fitGaussians(extractedBeadData);
			else
				means = fitFindPeaks(extractedBeadData);
			end
			
			% Filter means
			thresh = 8; % Not sure where this number comes from
			filteredMeans = filterMeans(means, counts, centers, thresh);
			
			numPeaks = size(filteredMeans, 1) - 1;	
			res_min = inf;
			opt = optimoptions('lsqcurvefit', 'Display', 'off');
			for hpi = 0 : (MAX_PEAKS - numPeaks)
				% Its a bit simpler to index peaks from the highest down,
				% since all beads have the largest populations, but some
				% don't have smaller ones. By flipping so that low indexes 
				% are high beads and vice-versa, we don't have to adjust 
				% the index for each bead population up or down.
				% --> hpi := highestPeakIndex
				
				highestPeak = BEAD_PEAKS - hpi;
				fprintf(1, 'Trying highest peak: %d\n', highestPeak);
				
				% Iterate over each channel and add up residuals
				% Reset rolling parameters each time for cleanliness
				res_all = zeros(1, numel(channels));
				fits = zeros(2, numel(channels));
				for chID = 1:numel(channels)
					
					% Extract vals 
					MEF = flipud(beadVals.(MEF_units{chID}));
					
					% Identify which bead pops to query
					beadPops = hpi + 1: (hpi + numPeaks);
					pointsMean = filteredMeans(2:end, chID); % Means already log10 transf
					pointsMEF = flipud(log10(MEF(beadPops)));
% 					pointsMEF = flipud((MEF(beadPops)));
					
					% Fit bead fluorescence to MEF values
					linearFunc = @(p, x) p(1) * x + p(2);
					if ismember('nonLinear', options)
						fit = lsqcurvefit(linearFunc, [1; 1], pointsMean, pointsMEF, [], [], opt);
					else
						fit = lsqcurvefit(@(p, x) linearFunc([1; p(2)], x), [1; 1], pointsMean, pointsMEF, [], [], opt);
					end
					fits(:, chID) = reshape(fit, [], 1);
					
					% Measure and collate residuals
					yResid = abs(pointsMEF - polyval(fit, pointsMean));
					ssResid = sum(yResid.^2);
					ssTotal = (length(pointsMEF) - 1) * var(pointsMEF);
					res_all(chID) = ssResid / ssTotal;
				end
				res = min(res_all);

				fprintf(1, 'Residuals: %.3f\n', res);
				
				% Determine if this is the best fit
				if ismember(sprintf('peak_%d', highestPeak), options)
					res = -inf;
					fprintf(1, 'Forcing Peak %d to be highest\n', highestPeak);
				end
				if (res < res_min)
					res_min = res;
					res_all_min = res_all;
					numPeaks_min = numPeaks;
					highestPeak_min = BEAD_PEAKS - hpi;
					beadPops_min = BEAD_PEAKS - beadPops + 1;
					fits_min = fits;
					means_min = 10.^filteredMeans(2:end, :);
% 						gmmMeans_min = means(2:end, :);
				end
				if (res_min == -inf) % The current peak was forced, skip the rest
					break
				end
			end
			
			% Extract finalized value
			rsq			= 1 - res_all_min;
			numPeaks	= numPeaks_min;
			hpi = highestPeak_min;
			beadPops	= beadPops_min;
			fits		= fits_min;
			meansLin	= means_min;
			meansBxp	= Transforms.lin2logicle(meansLin);
			
			fprintf(1, '\nFinished fitting! # Peaks = %d, Highest = Peak #%d\n', numPeaks, hpi);
			
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
					MEF = flipud(beadVals.(MEF_units{chID}));
					MEF = flipud(MEF(BEAD_PEAKS - beadPops + 1));
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
				end
			end
			
			% Convert fits to a table
			mefFits = array2table(fits);
			mefFits.Properties.RowNames = {'slope', 'zero'};
			mefFits.Properties.VariableNames = channels;
			
			
			% --- Helper Functions --- %
			
			
			function [beadData, beadVals] = zCheckInputs_calibrateMEF()
				
				validateattributes(beads, {'struct'}, {}, mfilename, 'beadsFilename', 1);
				minimumFields = {'filename', 'type', 'lot', 'date', 'cytometer'};
				missingFields  = setdiff(minimumFields, fieldnames(beads));
				assert(isempty(missingFields), 'Field not found: %s\n', missingFields{:})
				
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
				
				% Combine data from multiple bead samples if applicable 
				% (should just be a single file)
				if numel(beadData) > 1
					for ch = 1:numel(channels)
						channelData = [];
						for bdi = 1:numel(beadData)
							channelData = [channelData; beadData(bdi).(channels{ch}).raw]; %#ok<AGROW>
						end
						beadData(1).(channels{ch}).raw = channelData;
					end
					beadData(1).nObs = size(beadData.(channels{1}).raw, 1);
					beadData = beadData(1);
				end
								
				% Get bead MEF values
				beadVals = Transforms.getBeadVals(beads.type, beads.lot, channels);
			end
			
			
			function [extractedBeadData, figFits, counts, centers] = extractBeadData(beadData, showPlots)
				% Iterate over channels and extract bead data for requested channels
				% in order that they were given, find initial peaks for each channel
				% independently so we can figure out how many bead populations to fit.
				
				extractedBeadData = zeros(beadData.nObs, numel(channels));
				numBins = round(size(extractedBeadData, 1) / 100);
				figFits = struct();
				counts = struct();
				centers = struct();
				
				% Check CV of forward and side scatter to determine if gating is
				% needed to extract "good" bead population
				fscData = beadData.FSC_A.raw;
				sscData = beadData.SSC_A.raw;
				fscCV = std(fscData) / mean(fscData);
				sscCV = std(sscData) / mean(sscData);
				if true%(fscCV > 1.0 || sscCV > 2.0)
% 					warnString = 'High FSC/SSC Variance detected, please gate';
% 					warning(warnString)
					[gateIdx, ~, figFits.gate] = Gating.gatePolygon( ...
							fscData, sscData, struct('YScale', 'log', ...
									'XLabel', struct('String', 'FSC\_A'), ...
									'YLabel', struct('String', 'SSC\_A'), ...
									'Title', struct('String', 'Gate Beads Scatter'))); 
				else
					gateIdx = true(beadData.nObs, 1);
				end
				
				% Iterate over channels to extract data
				for ch = 1:numel(channels)
					
					% Extract channel data
					extractedBeadData(:, ch) = beadData.(channels{ch}).raw;
					
					% Create histogram w/ biexp data so that peaks are easier to find
					if showPlots
						figFits.(channels{ch}) = figure(); 
						subplot(2, 1, 1)
					end
					[counts.(channels{ch}), centers.(channels{ch})] = Plotting.biexhist( ...
								extractedBeadData(:, ch), numBins, showPlots);
				end
				
				% Only take positive data (that passes the gate if applicable)
				valid = (sum(extractedBeadData > 0, 2) == numel(channels));
				extractedBeadData = extractedBeadData((valid & gateIdx), :);
			end
			
			
			function [means, figManualPeaks] = fitGaussians(extractedBeadData)
				
				extrBeadDataBiex = Transforms.lin2logicle(extractedBeadData);
				
				% Have the user draw a set of lines where each point corresponds
				% with the center of a peak
				[~, sortIdx] = sort(max(extractedBeadData, [], 1));
				figManualPeaks = figure(); axG = gca();
				Plotting.densityplot(axG, ...
						extrBeadDataBiex(:, sortIdx(end)), ...
						extrBeadDataBiex(:, sortIdx(end-1)), ...
						5000, 'hist', ColorMap('parula'));
				title('Draw a shape connecting the peaks (in any order)')
				xlabel([strrep(channels{sortIdx(end)}, '_', '-'), ' (AFU)'])
				ylabel([strrep(channels{sortIdx(end - 1)}, '_', '-'), ' (AFU)'])
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
		
		
		function [meflFits, figFits] = calibrateMEFL(controlData, channels, colors, dataType, showPlots)
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
			%		colors			<cell, char> Name of fluorescent proteins and/or 
			%						fluorophores used in the experiment. 
			%						 - Should be equivalent in number to and match 
			%						   the order of 'channels'.
			%
			%		dataType		<char> The dataType to use for calculations
			%						 - This should be the post-MEF conversion 
			%						   compensated dataType (eg 'mComp')
			%
			%		showPlots		(optional) <logical> Flag to show fitting plots
			%
			%	Outputs
			%		meflFits		<struct> A struct array where each field is
			%						the channel name of a converted channel and the 
			%						value is the conversion factor to get MEFL units. 
			%       figFits         <struct> A struct of handles to all figures
            %						generated by the method. If 'showPlots' is not
			%						included as an optional input, this will return
			%						as empty.
			%
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			% 
			
			zCheckInputs_calibrateMEFL();
			FITC_IDX = find(strcmpi('FITC_A', channels));
			
			% Compute conversions
			meflFits = struct();
			figFits = struct();
			for chID = 1:numel(channels)
				if (chID == FITC_IDX) 
					% Skip FITC channel since it is the reference channel
					meflFits.(colors{chID}) = 1;
				else
					currChanMEF = Transforms.CHANNEL_MAP.(channels{chID});
					
					xdata = controlData(chID).(channels{chID}).(dataType);
					ydata = controlData(chID).FITC_A.(dataType);
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

					meflFits.(colors{chID}) = slope;
					
					if (exist('showPlots', 'var') && showPlots)
						
						figFits.(colors{chID}) = figure(); 
						ax = gca(); hold(ax, 'on');
						xrange = logspace(0, 9, 50);
						
						plot(ax, xdata, ydata, ...
							 '.', 'MarkerSize', 2)
						plot(ax, xrange, xrange * slope, ...
							 '-', 'LineWidth', 4)
						
						ax.YScale = 'log';
						ax.XScale = 'log';
						title(sprintf('Scale Factor: %.2f | R^2: %.3f', ...
								slope, rsq), 'fontsize', 14)
						ylabel([colors{FITC_IDX}, ' (MEFL)'])
						xlabel(sprintf('%s (%s)', strrep(colors{chID}, '_', '-'), currChanMEF));
					end
				end
			end
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_calibrateMEFL()
				validateattributes(controlData, {'struct'}, {}, mfilename, 'controlData', 1);
				validateattributes(channels, {'cell', 'char'}, {}, mfilename, 'channels', 2);
				if ischar(channels), channels = {channels}; end % Convert to cell for simplicity
				assert(numel(channels) == numel(controlData), ...
						'Number of channels and controls must be equal!\n');
				
				validateattributes(colors, {'cell', 'char'}, {}, mfilename, 'colors', 3);
				if ischar(colors), colors = {colors}; end % For simplicity
				assert(numel(colors) == numel(channels), ...
						'Number of colors and channels must be equal!\n');
				
				validChannels = fieldnames(controlData(1));
				badChannels = setdiff(channels, validChannels);
				assert(isempty(badChannels), 'Channel not in controlData: %s\n', badChannels{:});
				
				validatestring(dataType, fieldnames(controlData(1).(channels{1})), mfilename, 'dataType', 4);
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
			% Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
            
            % Process inputs
			channels = zCheckInputs_fcs2ABC(data, channelFits, dataType);
                        
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
			
            
			function channels = zCheckInputs_fcs2ABC(data, channelFits, dataType)
				
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
        
        
        function [abcFits, figFits] = calibrateABCfromMedians(beadsFnames, channel, hasBlank)
            % calibrateABCfromMedians creates antibody binding capacity (ABC) fits for an 
            % .fcs file of a particular bead sample supplied by the user using medians of 
            % each file. Note this is for Quantum Simply Cellular (R) beads, not rainbow 
            % beads (see calibrateMEF()). 
            %
            % Calibration plots and the table with ABC fits are saved in pwd().
            %
            %   Inputs: 
            %   
            %       beadFnames			A cell array of .fcs filenames which contain Ab-stained 
            %                           bead data. Can be given in any order (median is sorted).
            %
            %       channel             The channel name to calibrate bead data
            %
            %       hasBlank            (Optional) Boolean to indicate if the samples include the
            %                           blank population (ignore lowest peak if so). Default = false
            %
            %   Outputs:
            %
            %       abcFits				An Nx2 table of slopes/intercepts for linear conversion of 
            %                           N channels to ABC equivalents. Rows are labeled with 
            %                           channel names and columns are labeled with slope / zero
            %       
            %       figFits             A struct containing all figures generated by the function
            %
            % Written By
			% Ross Jones
			% jonesr18@mit.edu
            % Weiss Lab, MIT
			%
			% Update Log:
			%	2018-01-26:		Added figure handle output
            
            % Import data
			if ischar(beadsFnames), beadsFnames = {beadsFnames}; end
            beadData = FlowAnalysis.openFiles(beadsFnames{:});
            
            % Initialize as an array, convert to a table later
            abcFits = zeros(1, 2);
            
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
            figFits = figure('Position', [400 200 1000 400]);
            ax1 = subplot(1, 2, 1); hold(ax1, 'on');
            for i = 1:length(beadData)
                channelData = beadData(i).(channel).raw; 
				
				% Create histogram w/ biexp data so that peaks are easier to find
                [binCounts, binCenters] = Plotting.biexhist(channelData, 25, true);
                
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
            abcFits(1, :) = fit;
			
            fprintf(1, 'Finished fitting %s. # Peaks = %d\n', channel, length(peaks));
			
            % Convert fits to a table
            abcFits = array2table(abcFits);
            abcFits.Properties.VariableNames = {'slope', 'zero'};
            abcFits.Properties.RowNames = {channel};
			
            % Save channelFits
%             save([channel, '_ABC_Fit.mat'], 'channelFits')
        end
        
        
        function [abcFits, locs, peaks, figFits] = calibrateABC(beadsFilename, channel, largePeaks, hasBlank, locs, peaks)
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
            % Written By
			% Ross Jones
			% jonesr18@mit.edu
			% Weiss Lab, MIT
			% 
			% Update Log:
			%	2017-09-16:		Overhauled interface and processing
            
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
            abcFits = zeros(1, 2);

            % Extract channel data
            channelData = [];
            for i = 1:length(beadData)
                channelData = [channelData; beadData(i).(channel).raw]; %#ok<AGROW>
            end
            numBeads = length(channelData);
            
            % Plot figure to show fitting
            figFits = figure('Position', [400 200 1000 400]);
            ax1 = subplot(1, 2, 1); hold(ax1, 'on')
            numBins = round(numBeads / 100);
            [counts, centers] = Plotting.biexhist(channelData, numBins);
            
            % Find highest peaks
            if ~(exist('locs', 'var') && exist('peaks', 'var'))
                [peaks, locs] = findpeaks(counts);
                if (largePeaks)
                    maxpeak = max(peaks);
                    [peaks, locs] = findpeaks(counts, 'MinPeakHeight', maxpeak / 10);
                end
            end
            
            % Plot peaks on data
            plot(ax1, centers(locs), peaks, '*r')
            title(strrep(channel, '_', '-'))
            ylabel('COUNT')
            
            % Assume we have the P highest peaks and an unstained population
%             ABC = [68458, 281582, 621651];
            ABC = [5899, 68458, 281582, 621651];
			
            if (exist('hasBlank', 'var') && hasBlank)
                P = min(length(locs) - 1, numel(ABC));
            else
                P = min(length(locs), numel(ABC));
            end
            ABC = ABC(end-P+1:end);
            peakIntensity = Transforms.logicle2lin(centers(locs));
			peakIntensity = peakIntensity(end-P+1:end);
            
            % Calculate fit
            fit = polyfit(peakIntensity, log10(ABC), 1);
            
            ax2 = subplot(1, 2, 2); hold(ax2, 'on')
            semilogy(ax2, peakIntensity, ABC, 'ob')
            semilogy(ax2, peakIntensity, 10.^polyval(fit, peakRelChannel), 'r-')
            
            % Calculate R^2 value
            yresid = ABC - 10.^polyval(fit, peakIntensity);
            SSresid = sum(yresid.^2);
            SStotal = (length(ABC) - 1) * var(ABC);
            rsq = 1 - SSresid / SStotal;
            title(sprintf('R^2 = %.3f', rsq))
            axis([0 2^18 1 1e6])
            xlabel('REL CH #')
            ylabel('ABC')
            hold off

            % Assign fit data
            abcFits(1, :) = fit;

            fprintf(1, 'Finished fitting %s. # Peaks = %d\n', channel, length(peaks));
            savefig(figPeaks, [channel '_ABC_Fit.fig'])

            % Convert fits to a table
            abcFits = array2table(abcFits);
            abcFits.Properties.VariableNames = {'slope', 'zero'};
            abcFits.Properties.RowNames = {channel};

            % Save channelFits
            save([channel, '_ABC_Fit.mat'], 'channelFits')
		end
        
		
		function params = checkLogicleParams(doMEF, params)
			% Simple function for checking and setting default logicle parameters
			
			% Set parameters
			if ~isfield(params, 'T'), params.T = 2^18;	end
			if ~isfield(params, 'M'), params.M = 4.5;	end
			if ~isfield(params, 'r'), params.r = -150;	end
			if doMEF 
				if ~isfield(params, 'MEF')
					params.MEF = Transforms.MEF_CONVERSION_FACTOR;
				end
			else
				% lin2logicle and logicle2lin use params.MEF for scaling, so if
				% we are not doing a MEF conversion, set to 1. 
				params.MEF = 1;
			end
		end
		
		
		function params = checkLogicleParams2(params)
			% Simple function for checking and setting default logicle parameters
			
			% Set parameters
			if ~isfield(params, 'T'), params.T = 2^18;	end
			if ~isfield(params, 'M'), params.M = 4.5;	end
			if ~isfield(params, 'r'), params.r = -150;	end
		end

		
    end
    
end