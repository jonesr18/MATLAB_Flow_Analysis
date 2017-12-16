classdef FlowData < handle
	% A data structure for managing flow cytometry data. Uses and references
	% several other classes and functions in the repository. 
	%
	%	Visible Properties
	%
	%		name			<char>		Experiment name
	%		date			<char>		Experiment start date
	%
	%		numSamples		 <numeric>	The number of data samples
	%		numCells		 <array>	The number of cells in each sample
	%		sampleData		 <struct>	Sample fluorescence and gate data in standard struct
	%		sampleMap		 <table>	Experimental information for samples
	%		dataTypes		 <cell>		Cell array of data types 
	%									('raw', 'mComp', 'mefl', etc)
	%		gateNames		 <cell>		Cell array of gate names (strings)
	%		gatePolygons	 <struct>	Mapping between gate names and polygons for sampleData
	%		channels		 <cell>		Cell array of channel names
	%		controlData		 <struct>	Similar to sampleData but for controls.
	%									The order of controls should be the same
	%									as their corresponding channels. 
	%		fitParams		 <struct>	Struct containing matrix compensation fitting results.
	%		compDataType	 <char>		Data type used for compensation
	%		beadFitsControls <table>	Table containing MEF unit fits for controls
	%		meflConversions  <struct>	Mapping between channel names to MEF-MEFL conversion factors 
	%		beadFitsSamples	 <table>	Table containing MEF unit fits for samples
	%		bins			 <cell>		Cell array where each element corresponds
	%									with a data sample and contains an N-dim
	%									M-element/dim cell array with an array of
	%									numerical cell indexes (IDs) in each bin. 
	%										N = numel(binChannels)
	%										M = numBins
	%										IDs will depend on binDataType & binGate
	%		binInputs		 <struct>	Struct w/ bin channels as fields and edges as values
	%		binDataType		 <char>		Binned dataType ('raw', 'mComp', 'mefl', etc)
	%		binGate			 <char>		Gate name for cell population to bin
	%		binStats		 <struct>	Struct w/ most recently generated bin statistics.
	%									Each field is a different measurement (median, 
	%									stdev, etc) and each value is an N-dimensional 
	%									matrix containing the measured values. The
	%									dimensionality is determined by which channels
	%									are requested in self.getBinStats()
	%
	%	Public Methods
	%
	%		FlowData
	%		addControls
	%		gate
	%		convertToMEF
	%		compensate
	%		convertToMEFL
	%		bin
	%		getBinStats
	%		getValues
	%		getSampleIDs
	%		slice
	%		threshGate
	%		crossGates
	
	%#ok<*AGROW>	
	
	properties (SetObservable, Hidden)
		test = '';
	end
	
	
	properties (SetAccess = private)
		name = '';					% Experiment name
		date = '';					% Experiment start date
		cytometer = '';				% Experiment cytometer used
		
		numSamples = 0;				% The number of data samples
		numCells = [];				% The number of cells in each sample
		sampleData = struct();		% Sample fluorescence and gate data in standard struct
		sampleMap = table();		% Experimental information for samples
		
		gateNames = {};				% Cell array of gate names (strings)
		gatePolygons = struct();	% Mapping between gate names and polygons for sampleData
		
		dataTypes = {'raw'};		% Cell array of data types ('raw', 'mComp', 'mefl', etc)
		channels = {};				% Cell array of channel names
		controlData = struct();		% Similar to sampleData but for controls
		fitParams = struct();		% Struct containing matrix compensation fitting results
		compDataType = '';			% Data type used for compensation
		
		beadFitsControls = table();	% Table containing MEF unit fits for controls
		beadFitsSamples = table();	% Table containing MEF unit fits for samples
		meflConversions = struct(); % Mapping between channel names to MEF-MEFL conversion factors 
		
		bins = {};					% Cell array where each element corresponds with a data sample
		binInputs = struct();		% Struct w/ bin channels as fields and edges as values
		binDataType = '';			% Binned dataType ('raw', 'mComp', 'mefl', etc)
		binGate = '';				% Gate name for cell population to bin
		binStats = struct();		% Most recently generated bin statistics
	end
	
	
	properties (Access = private)
		sampleDataScatter = struct();	% Contains sampleData forward/side scatter values
		controlDataScatter = struct();	% Contains controlData forward/side scatter values
		onlyP1 = false;					% Flag to only gate w/ P1 (see self.gate() )
		
		controlsAdded = false;			% Tracks if controls have been added
		gated = false;					% Tracks if gating has been performed
		mefConverted = false;			% Tracks if mef conversions have been performed
		compensated = false;			% Tracks if compensation has been performed
		meflConverted = false;			% Tracks if mefl conversion has been performed
		binned = false;					% Tracks if binning has been performed
	end
	
	
	properties (Access = private, Constant)		
		SHORT_COLORS = struct( ...
				'BUV_396_A', 'V', ...
				'Pacific_Blue_A', 'B', ...
				'FITC_A', 'Y', ...
				'PE_A', 'O', ...
				'PE_Texas_Red_A', 'R', ...
				'PE_TxRed_YG_A', 'R', ...
				'APC_Cy7_A', 'I');
	end
	
	
	methods (Access = public)
		
		function self = FlowData(dataFnames, channels, exptDetails, sampleMapFname)
			% Initializes the FlowData object 
			% by importing data from the given files, which should correspond
			% with the given sample map and contain data in the given channels.
			%
			%	self = FlowData(dataFilenames, channels, sampleMap)
			%
			%	Inputs
			%		dataFnames		<cell, char> Data filenames to import
			%						** The filenames will be sorted on their
			%						last 3 characters (excluding '.fcs'). 
			%						The sorted order should coincides with 
			%						sample numbers in sampleMap.
			%
			%		channels		<cell, char> The color channel(s) corresponding
			%						with the desired subset of data in dataStruct.
			%						**Must be a field of dataStruct
			%						**FSC/SSC are automatically taken
			%
			%		exptDetails		<struct> A struct with the following fields
			%						recording experimental details:
			%							name  |  <char>
			%							date  |  <char>
			%
			%		sampleMapFname	<char> The name of the .txt file containing 
			%						treatments information for each sample. 
			%						**See example file in source folder.
			%
			%	Outputs
			%		self			A handle to the object
			
			% When loading, we need to make a blank version of the object
			if strcmpi(dataFnames, 'load'), return, end
			
			dataStruct = checkInputs();
			
			% Extract experiment details
			self.date = exptDetails.date;
			self.name = exptDetails.name;
			self.cytometer = exptDetails.cytometer;
			
			% Extract sample map
			sampleMap = readtable(sampleMapFname, 'delimiter', '\t');
			assert(height(sampleMap) == numel(dataStruct), ...
				'Sample map (%d) contains incorrect number of rows (%d)', ...
				height(sampleMap), numel(dataStruct))
			if ismember('Replicate', sampleMap.Properties.VariableNames)
				numReplicates = max(sampleMap.Replicate);
				ids = zeros(numReplicates, height(sampleMap) / numReplicates);
				for r = 1:numReplicates
					ids(r, :) = find(sampleMap.Replicate == r);
				end
				% Wean sampleMap to only be have one replicate
				sampleMap = sampleMap(sampleMap.Replicate == 1, :);
			else
				numReplicates = 1;
				ids = 1:height(sampleMap);
			end
			self.sampleMap = sampleMap;
			
			% Extract data from the given channels
			self.numSamples = height(sampleMap);
			self.numCells = zeros(1, self.numSamples);
			iii = 0;
			for id = ids
				
				iii = iii + 1;
				nObs = 0;
				for r = 1:numReplicates
					nObs = nObs + dataStruct(id(r)).nObs;
				end
				self.numCells(iii) = nObs;
				
				% Extract desired color channels
				for ch = channels
					sd = [];
					for r = 1:numReplicates
						sd = [sd; dataStruct(id(r)).(ch{:}).raw];
					end
					self.sampleData(iii).(ch{:}).raw = sd;
				end
				self.sampleData(iii).nObs = nObs;
				
				% Extract scatter channels
				for ch = Gating.SCATTER_CHANNELS
					sds = [];
					for r = 1:numReplicates
						sds = [sds; dataStruct(id(r)).(ch{:}).raw];
					end
					self.sampleDataScatter(iii).(ch{:}).raw = sds;
				end
				self.sampleDataScatter(iii).nObs = nObs;
			end
			self.channels = channels;
			
			% Add listeners for settable public properties
% 			addlistener(self, 'test', 'PostSet', @self.handlePropEvents);
% 			addlistener(self, {'binInputs', 'binDataType', 'binGate'}, ...
% 				'PostSet', @self.handlePropEvents);
			
			fprintf(1, 'Finished constructing FlowData object\n')
			
			
			% -- Helper functions -- %
			
			
			function dataStruct = checkInputs()
				validateattributes(dataFnames, {'cell', 'char'}, {}, mfilename, 'dataFilenames', 1);
				validateattributes(channels, {'cell', 'char'}, {}, mfilename, 'channels', 2);
				validateattributes(exptDetails, {'struct'}, {}, mfilename, 'exptDetails', 3);
				validateattributes(sampleMapFname, {'char'}, {}, mfilename, 'sampleMap', 4);
				
				% Convert channels to cell array if single char value is given
				if ischar(channels), channels = {channels}; end
				channels = reshape(channels, 1, []); % Ensure row vector
				
				% Check required experiment details are present
				requiredFields = {'date', 'name', 'cytometer'};
				missingFields = setdiff(requiredFields, fieldnames(exptDetails));
				assert(isempty(missingFields), 'Experiment details missing field: %s\n', missingFields{:});
				
				% Check sampleMap is a real file
				assert(logical(exist(sampleMapFname, 'file')), ...
					'File not found: %s\n', sampleMapFname)
				
				% Import data
				dataStruct = FlowAnalysis.openFiles(dataFnames{:});
				
				% Check channels are present in dataStruct
				badChannels = setdiff(channels, fieldnames(dataStruct(1)));
				assert(isempty(badChannels), ...
					'Channel not in dataStruct: %s\n', badChannels{:});
			end
		end
		
		
		function addControls(self, wildTypeFilename, singleColorFilenames, twoColorFilenames)
			% Adds wild-type, single-color, and two-color (optional) data to the dataset
			% so that we can do compensation (single-colors) and MEFL conversion (two-colors).
			%
			% The method generates self.controlData, a struct array where single-color  
			% data from channel X is in position X, two-color data from channel X is
			% in position 2*X (if applicable), and wild-type data is in the last
			% position (regardless of the presence of two-color data). 
			% 
			%	self.addControls(wildTypeFilename, singleColorFilenames, twoColorFilenames)
			%
			%	Inputs
			%		wildTypeFilenames		<cell, char> Wild-type cell data
			%		
			%		singleColorFilenames	<cell, char> Single-color controls data files
			%								** Order of colors should coincide with
			%								the order of FlowData.channels
			%
			%		twoColorFilenames		(optional) <cell, char> Two-color controls data files
			%								** As with single colors, the order should match
			%								those in self.channels, but with no yellow/green 
			%								file (since all other colors are converted to MEFL 
			%								units using the FITC channel).
			
			[wildTypeData, singleColorData, twoColorData] = checkInputs_addControls(self);
			FITC_IDX = find(strcmpi('FITC_A', self.channels));
			
			% Extract data
			self.controlData = extractData([self.channels, {'nObs'}], ...
						wildTypeData, singleColorData, twoColorData, FITC_IDX);
			
			% Extract scatter data
			self.controlDataScatter = extractData([Gating.SCATTER_CHANNELS, {'nObs'}], ...
						wildTypeData, singleColorData, twoColorData, FITC_IDX);
			
			self.controlsAdded = true;
			fprintf(1, 'Finished adding controls\n');
			
			
			% --- Helper Functions --- %
			
			
			function [wildTypeData, singleColorData, twoColorData] = checkInputs_addControls(self)
				
				validateattributes(wildTypeFilename, {'cell', 'char'}, {}, mfilename, 'wildTypeFilename', 1);
				validateattributes(singleColorFilenames, {'cell', 'char'}, {}, mfilename, 'singleColorFilenames', 2);

				% Convert to cell arrays if necessary for convenience
				if ischar(wildTypeFilename), wildTypeFilename = {wildTypeFilename}; end
				if ischar(singleColorFilenames), singleColorFilenames = {singleColorFilenames}; end

				% Check number of scFiles
				assert(numel(singleColorFilenames) == numel(self.channels), ...
					'Incorrect number of single color controls');

				% Open files
				wildTypeData = FlowAnalysis.openFiles(wildTypeFilename{:});
				singleColorData = FlowAnalysis.openFiles(singleColorFilenames{:});

				% Add twoColorData if applicable
				if exist('twoColorFilenames', 'var')
					validateattributes(twoColorFilenames, {'cell', 'char'}, {}, mfilename, 'twoColorFilenames', 3);
					if ischar(twoColorFilenames), twoColorFilenames = {twoColorFilenames}; end
					assert(numel(twoColorFilenames) == sum(~strcmpi('FITC_A', self.channels)), ...
						'Incorrect number of two color controls');
					twoColorData = FlowAnalysis.openFiles(twoColorFilenames{:});
				else
					twoColorData = [];
				end
			end
			
			
			function outData = extractData(channels, wtData, scData, tcData, FITC_IDX)
				% Generates the new controlData struct based on the given
				% individual structs and the desired channels
				%
				%	The FITC_IDX input is needed to tell the function which
				%	channel ID to skip when extracting two-color data. There is
				%	no two-color data for FITC since the FITC channel itself is
				%	the reference for the other fluorescent proteins. 
				
				for ch = channels
					for sc = 1:numel(scData)
						outData(sc).(ch{:}) = scData(sc).(ch{:});
					end
					
					tc = 0;
					if ~isempty(tcData)
						tcIdx = 0;
						for tc = 1:numel(scData) % Go over length of scData since they should match
							if tc == FITC_IDX
								% No two color controls for FITC channel, since it
								% is the reference color for MEFL conversion
								outData(sc + tc).(ch{:}) = [];
							else
								tcIdx = tcIdx + 1; % For indexing tcData separate of tc iterator
								outData(sc + tc).(ch{:}) = tcData(tcIdx).(ch{:});
							end
						end
					end
					outData(sc + tc + 1).(ch{:}) = wtData.(ch{:});
				end
			end
		end
		
		
		function gate(self, onlyP1)
			% Creates gates for the data using standard gating (see Gating.m)
			%
			%	self.gate(onlyP1)
			%
			%	onlyP1 = true adds new gates: {'P1'}
			%	onlyP2 = false adds new gates: {'P1', 'P2', 'P3'}
			%
			%	Inputs
			%		onlyP1		<logical> Logical flag to only do P1 gating
			%						(FSC_A vs SSC_A)
			
			assert(self.controlsAdded, 'Controls must be added before gating!\n');
			
			% Determine gating options
			onlyP1 = (~exist('onlyP1', 'var') || all(logical(onlyP1)));
			self.onlyP1 = onlyP1;
			if strcmpi(self.cytometer, 'Koch-LSRII-HTS2'), swap = true; else, swap = false; end
			
			% Setup new directory for gates
			gateDir = ['Gating', filesep];
			if ~exist(gateDir, 'file')
				mkdir(gateDir)
			end
			
			% Check if gates have already been made for this data
			% --> If so, load them and immediately return!
			gatesSaveName = [self.date, '_', self.name];
			gatesFname = [gateDir, gatesSaveName '_GatePolygons.mat'];
			gatesFnameControls = [gateDir, 'Controls', '_GatePolygons.mat'];
			if exist(gatesFname, 'file') 
				% Load existing sample gates
				load(gatesFname);
			else
				% Do manual sample gating
				[gateP1s, gateP2s, gateP3s] = Gating.standardGating(self.sampleDataScatter, onlyP1, swap);
				save(gatesFname, 'gateP1s', 'gateP2s', 'gateP3s');
			end
			if exist(gatesFnameControls, 'file')
				% Load existing control gates
				load(gatesFnameControls);
			else
				% Do manual control gating
				[gateP1c, gateP2c, gateP3c]	= Gating.standardGating(self.controlDataScatter, onlyP1, swap);
				save(gatesFnameControls, 'gateP1c', 'gateP2c', 'gateP3c');
			end
			
			% Only extract polygons and such for sampleData, not controlData
			self.gatePolygons.P1 = gateP1s;
			self.addGates('P1');
			if ~onlyP1
				self.gatePolygons.P2 = gateP2s;
				self.gatePolygons.P3 = gateP3s;
				self.addGates({'P2', 'P3'});
			end
			
			% Apply gate polygons to scatter data 
			self.controlDataScatter = Gating.applyStandardGates(self.controlDataScatter, gateP1c, gateP2c, gateP3c, swap);
			self.sampleDataScatter	= Gating.applyStandardGates(self.sampleDataScatter,  gateP1s, gateP2s, gateP3s, swap);
			
			% Transfer gate logicals to extrenally accesible data
			for cd = 1:numel(self.controlData)
				self.controlData(cd).gates = self.controlDataScatter(cd).gates;
			end
			for sd = 1:self.numSamples
				self.sampleData(sd).gates = self.sampleDataScatter(sd).gates;
			end
			
			fprintf(1, 'Finished standard gating\n')
		end
		
		
		function convertToMEF(self, beadsControls, beadsSamples, options)
			% Calibrate the data to standardized bead-based MEF units.
			% The controls must already have been added (see addControls())
			%
			%	self.convertToMEF(beadsFilename)
			%	
			%	Default settings force a linear fit (in log space), ie the
			%	log-log fit is forces to have a slope of 1, which ensures that
			%	the MEF conversion is a simple scaling function of fluorescent
			%	to MEF units. 
			%
			%	Adds new dataTypes: {'mef'}
			%	Adds new gates: {'nneg', 'P3_nneg'}
			%
			%	Inputs
			%		beadsControls	<struct> A struct with the following fields:
			%			filename		The .fcs file containing bead data
			%							corresponding with this experiment
			%			type			The name of the type of bead (eg 'RCP-30-5A')
			%			lot				The bead production lot (eg 'AH01')
			%		beadsSamples	<struct> A struct with the same fields as above,
			%						but for the samples rather than controls
			%		options			<cell> Contains optional string flags:
			%			showPlots		Flag to show fitted plots
			%			nonLinear		Flag to *not* force linear bead fits
			%			fitGaussian		Flag to use multi-dimensional gaussian
			%							fitting rather than the single-dim 
			%							peak-finding algorithm
			%			peak_X			Forces peak X (where X is replaced with
			%							a number) to be the highest bead peak,
			%							rather than having the function iterate
			%							through all possible high peaks automaticallly
			%			noFilter		Do not filter out peaks by height

			
			assert(self.controlsAdded, 'Controls must be added before converting to MEF units!\n');
			
			% Check inputs
			checkInputs_convertToMEF();
			
			% Get MEF fits 
			fitsControls = Transforms.calibrateMEF(beadsControls, self.channels, options);
			if isequaln(beadsControls, beadsSamples)
				% If the bead properties are the same, then skip the fitting for
				% samples' beads and just use the controls' beads. This is for 
				% the case where samples and controls are run the same day.
				fitsSamples = fitsControls;
			else
				fitsSamples = Transforms.calibrateMEF(beadsSamples, self.channels, options);
			end
			
			% Save fit information
			self.beadFitsControls = fitsControls;
			self.beadFitsSamples = fitsSamples;
			
			% Apply calibration to controls
			self.controlData = Transforms.fcs2MEF(self.controlData, fitsControls, 'raw');
			
			% Apply calibration to data
			self.sampleData = Transforms.fcs2MEF(self.sampleData, fitsSamples, 'raw');
			
			self.addDataTypes('mef');
			self.addGates('nneg');
			if self.onlyP1
				self.crossGates({'P1', 'nneg'}, 'and');
			else
				self.crossGates({'P3', 'nneg'}, 'and');
			end
			
			self.mefConverted = true;
			fprintf(1, 'Finished converting to MEF\n')
			
			
			% --- Helper Functions --- %
			
			
			function checkInputs_convertToMEF()
				
				validateattributes(beadsControls, {'struct'}, {}, mfilename, 'beadsFilename', 1);
				validateattributes(beadsSamples, {'struct'}, {}, mfilename, 'beadsFilename', 2);
				
				% No need to check bead fields, since that is taken care of by
				% Transforms.calibrateMEF()
				assert(logical(exist(beadsControls.filename, 'file')), ...
					'File does not exist: %s\n', beadsControls.filename)
				assert(logical(exist(beadsSamples.filename, 'file')), ...
					'File does not exist: %s\n', beadsSamples.filename)
				
				assert(strcmpi(beadsControls.cytometer, self.cytometer), ...
					'Control beads cytometer does not match experiment cytometer!')
				assert(strcmpi(beadsSamples.cytometer, self.cytometer), ...
					'Sample beads cytometer does not match experiment cytometer!')
				
				if ~exist('options', 'var')
					options = {};
				end
			end
		end
		
		
		function compensate(self, method, dataType, gate, plotsOn)
			% Applies compensation to the data
			%
			%	self.compensate(method, plotsOn)
			%
			%	method = 'scComp' adds new dataTypes: {'scComp'}
			%	method = 'mComp' adds new dataTypes: {'mComp', 'afs'}
			%
			%	Inputs
			%		method			<char> Indicates which compensation routine to use. 
			%							'scComp'	piecewise linear
			%							'mComp'		matrix-based (incl autorfluor subtract)
			%
			%		dataType		<char> Indicates which data type to use.
			%						Can be any dataType in self.dataTypes, but
			%						preferrably 'mef'.
			%
			%		gate			<char> Indicates which gate to use for.
			%						Can be any gate in self.gateNames, but
			%						preferrably 'P1', or 'P3' if onlyP1 = false
			%
			%		plotsOn			(optional) <logical> Set to TRUE to show the 
			%						compensation function's generated plots. 
			
			% Check inputs
			validatestring(method, {'scComp', 'mComp'}, mfilename, 'method', 1);
			validatestring(dataType, self.dataTypes, mfilename, 'dataType', 2);
			validatestring(gate, self.gateNames, mfilename, 'gate', 3);
			plotsOn = (exist('plotsOn', 'var') && all(logical(plotsOn)));
			
			switch method
				case 'scComp'
					self.sampleData = Compensation.compensateBatchSC( ...
						self.controlData(end), self.controlData(1:numel(self.channels)), ...
						self.channels, self.sampleData, self.channels, ...
						dataType, gate, plotsOn);
					self.addDataTypes('scComp');
				
				case 'mComp'
					[self.sampleData, self.controlData, self.fitParams] = ...
						Compensation.compensateMatrixBatch( ...
							self.sampleData, self.controlData, self.channels, ...
							dataType, gate, plotsOn);
					self.addDataTypes({'afs', 'mComp'});
			end
			
			self.compensated = true;
			self.compDataType = dataType;
			fprintf(1, 'Finished compensation\n')
		end
		
		
		function convertToMEFL(self, showPlots)
			% Converts each channel to MEFL units using the compensated MEF units. 
			% The two-color controls are utilized to get ratios between each MEF
			% unit and MEFLs. Conversion factors are stored in self.meflConversions.
			%
			%	self.convertToMEFL(showPlots);
			%
			%	Adds new dataTypes: {'mefl'}
			%
			%	Inputs
			%		showPlots	(optional) <logical> Flag to show conversion plots
			
			% Check pre-requisites for running
			assert(self.controlsAdded, 'Controls must be added before converting to MEFL units!\n');
			assert(self.mefConverted, 'MEF conversion must be run before converting to MEFL units!\n');
			assert(self.compensated, 'Compensation must be run before converting to MEFL units!\n');
			assert(strcmpi(self.compDataType, 'mef'), 'Compensation must be run on MEF-converted data!\n')
			
			% Choose dataType to convert
			if ismember('mComp', self.dataTypes)
				dataType = 'mComp';
			elseif ismember('scComp', self.dataTypes)
				dataType = 'scComp';
			else
				error('Compensated data not found!\n')
			end
			
			% Compute mefl conversions
			meflSaveName = self.name;
			tcData = self.controlData(numel(self.channels) + 1 : 2 * numel(self.channels));
			meflFits = Transforms.calibrateMEFL(tcData, self.channels, dataType, meflSaveName, showPlots);
			
			% Add converted data as 'mefl' data type
			for ch = self.channels
				% Add MEFLs for controls
				for i = 1:numel(self.controlData)
					if isempty(self.controlData(i).(ch{:})), continue, end % Some tcData will be empty 
					self.controlData(i).(ch{:}).mefl = self.controlData(i).(ch{:}).(dataType) * meflFits.(ch{:});
				end
				
				% Add MEFLs for sample
				for i = 1:self.numSamples
					self.sampleData(i).(ch{:}).mefl = self.sampleData(i).(ch{:}).(dataType) * meflFits.(ch{:});
				end
			end
			
			self.meflConversions = meflFits; 
			self.addDataTypes('mefl');
			self.meflConverted = true;
			fprintf(1, 'Finished converting to MEFL\n');
		end
		
		
		function binStats = bin(self, binInputs, binDataType, binGate)
			% Bins the data by the input channels into the given number of bins 
			%
			%	self.bin(binEdges, binChannels, binDataType, binGate)
			%
			%	Inputs
			%		binInputs		<struct> A struct with channel names as keys and 
			%						bin edges as values. The channel names must match 
			%						a field in the data struct with the subfield 'raw'. 
			%						The struct tells the function which channels to bin 
			%						on and where to draw the bins in each dimension. 
			%		
			%		binDataType		<char> The cell dataType to use (eg 'mefl', 'mComp')
			%
			%		binGate			(optional) <char> The gated cell population 
			%						to use for binning (default: 'P3')
			%
			%	Outputs
			%		binStats		<struct> A nested struct where each element corresponds
			%						with a sample and has the following heirarchy: 
			%							binStats --> channel --> statistic --> calculation
			%							Statistics:
			%								numCells, 10th, 50th (median), 90th %iles, 
			%								mean, geomean, stdev, geostdev, sem, semb*
			%								 * semb = bootstrapped SEM, currently
			%								   <not done because it is too slow>
			%
			%	Implementation notes:
			%		The binning method (FlowAnalysis.simpleBin()) operates in 
			%		roughly O(B*N) time where B is the number of bins and N is 
			%		the total number of cells. The "N" term is "fast" though, 
			%		since we take advantage of MATLAB's optimized matrix methods. 
			%
			%		An alternative method operates in roughly O(log(B)*N) time,
			%		but doesn't use the optimized matrix methods and is thus
			%		actually slower. 
			%
			%		The reason to do this maximally slow binning method is so
			%		that accessing which cells are in each bin operates in O(1)
			%		time, as we can directly query a bin to see which cells are
			%		inside of it, making all operations _after_ binning much 
			%		faster, overall improving the speed of analysis. 
			
			[binChannels, binEdges] = checkInputs_bin(self);
			
			for i = 1:self.numSamples
				
				sliceParams = struct( ...
					'channels', {binChannels}, ...
					'dataType', binDataType);
				
				slicedData = self.slice(i, sliceParams);
				
				% Transform to logicle space based on the dataType
				if (strcmpi(binDataType, 'raw') || ...
						(ismember(binDataType, {'afs', 'scComp', 'mComp'}) && ...
						 strcmpi(self.compDataType, 'raw')))
					% Raw data (comp or not)
					slicedData = Transforms.lin2logicle(slicedData);
% 					disp('Data tranformed to log')
				else
					% MEF or MEFL transformed data (comp or not)
					slicedData = Transforms.lin2logicleMEF(slicedData);
% 					disp('Data transformed to logMEF')
				end
				
				self.bins{i} = FlowAnalysis.simpleBin(slicedData, binEdges);
				
			end
			
			binStats = self.computeBinStats(binDataType, binGate);
			
			self.binInputs = binInputs;
			self.binDataType = binDataType;
			self.binGate = binGate;
			self.binned = true;
			self.binStats = binStats;
			fprintf(1, 'Finished binning\n')
			
			
			% --- Helper Functions --- %
			
			
			function [binChannels, binEdges] = checkInputs_bin(self)
				% Validates that the given bin properties are ok, then sets the
				% object's properties themselves if all are ok.

				% Check properties
				validateattributes(binInputs, {'struct'}, {}, mfilename, 'inputs', 1);
				binChannels = reshape(fieldnames(binInputs), 1, []);
				validateattributes(binChannels, {'cell', 'char'}, {}, mfilename, 'binChannels');
				badChannels = setdiff(binChannels, self.channels);
				assert(isempty(badChannels), ...
						'Channel not allowed: %s\n', badChannels{:});
				
				binEdges = cell(1, numel(binChannels));
				for bc = 1:numel(binChannels)
					assert(numel(binInputs.(binChannels{bc})) > 1, ...
							'Must have more than one bin edge to define a bin!');
					binEdges{bc} = binInputs.(binChannels{bc});
				end
				
				validateattributes(binDataType, {'char'}, {}, mfilename, 'binDataType', 2);
				assert(any(strcmp(binDataType, self.dataTypes)), ...
						'Bin data type does not match any existing data types: %s\n', binDataType);

				validateattributes(binGate, {'char'}, {}, mfilename, 'binGate', 3);
				assert(any(strcmp(binGate, self.gateNames)), ...
						'Gate does not exist in data: %s\n', binGate);
			end
		end
		
		
		function binStats = computeBinStats(self, dataType, gate)
			% Computes bin statistics based on the given bin dataType and bin gate. 
			%
			%	Inputs
			%		dataType	<char> Data type to compute stats with
			%		gate		<char> Gate to cross with the bins
			%
			%	Outputs
			%		binStats	<struct> A nested struct where each element corresponds
			%					with a sample and has the following heirarchy: 
			%						binStats --> channel --> statistic --> calculation
			%						Statistics:
			%							numCells, 10th, 50th (median), 90th %iles, 
			%							mean, geomean, stdev, geostdev, sem, semb*
			%							 * semb = bootstrapped SEM, currently
			%							   <not done because it is too slow>
			
			checkInputs_computeBinStats(self);
			
			statTypes = {'p10', 'p50', 'p90', 'mean', 'geomean', 'stdev', 'geostdev', 'sem'};
			
			for i = 1:self.numSamples
				
				cellsIdx = 1:self.numCells(i);
				passGate = cellsIdx(self.sampleData(i).gates.(gate));
				binStats(i).numCells = zeros(size(self.bins{i}));
				
				% Setup stats vectors
				for ch = 1:numel(self.channels)
					
					channel = self.channels{ch};
					
					% Set up binStats struct
					for st = statTypes
						binStats(i).(channel).(st{:}) = zeros(size(self.bins{i}));
					end
				end
				
				for b = 1:numel(self.bins{i})
					
					% Extract data
					inBin = intersect(self.bins{i}{b}, passGate);
					binStats(i).numCells(b) = numel(inBin);
					
					for ch = 1:numel(self.channels)
					
						channel = self.channels{ch};
						
						% Get data
						dataInBin = self.sampleData(i).(channel).(dataType)(inBin);
						dataInBin = dataInBin(~isnan(dataInBin));
						posData = (dataInBin > 0);
						
						% Calculate
						prctiles = prctile(dataInBin, [10, 50, 90]);
						binStats(i).(channel).p10(b) = prctiles(1);
						binStats(i).(channel).p50(b) = prctiles(2);
						binStats(i).(channel).p90(b) = prctiles(3);
						binStats(i).(channel).mean(b) = mean(dataInBin);
						binStats(i).(channel).geomean(b) = geomean(dataInBin(posData));
						binStats(i).(channel).stdev(b) = std(dataInBin);
						binStats(i).(channel).geostdev(b) = geostd(dataInBin(posData));
						binStats(i).(channel).sem(b) = std(dataInBin) / sqrt(numel(dataInBin));
% 						binStats(i).(channel).semb(b) = semBootstrap(dataInBin); % SUUUUPER SLOW
% 						binStats(i).(channel).CI95 = ci95(dataInBin);
					end
				end
			end
			
			
			% --- Helper Functions --- %
			
			
			function checkInputs_computeBinStats(self)
				
				validateattributes(dataType, {'char'}, {}, mfilename, 'binDataType', 1);
				assert(any(strcmp(dataType, self.dataTypes)), ...
						'Bin data type does not match any existing data types: %s\n', dataType);

				validateattributes(gate, {'char'}, {}, mfilename, 'binGate', 2);
				assert(any(strcmp(gate, self.gateNames)), ...
						'Gate does not exist in data: %s\n', gate);
				
			end
		end
		
		
		function values = getValues(self, varargin)
			% Returns unique values for each given experimental parameter.
			%
			%	values = self.getValues(paramters)
			%
			%	Inputs
			%		varargin	<char, cell> A cell array of strings or individual 
			%					string inputs of parameter names corresponding with 
			%					paramters in sampleMap. 
			%	
			%	Outputs
			%		values		<struct> a struct where each parameter is a field
			%					containing each unique parameter value as found
			%					in sampleMap. Values are produced as row vectors.
			
			% Check parameters
			parameters = {};
			for i = 1:numel(varargin)
				
				% Check each input to see if it is a single char or a cell of
				% filenames, since the method should take in variable inputs
				if iscell(varargin{i})
					parameters = [parameters, varargin{i}];
				else
					validateattributes(varargin{i}, {'char'}, {}, mfilename, 'parameters', 1);
					parameters = [parameters, varargin(i)];
				end
			end
			badParameters = setdiff(parameters, self.sampleMap.Properties.VariableNames);
			assert(isempty(badParameters), 'Parameter not recognized: %s\n', badParameters{:});
			
			% Extract unique values for each parameter, force to be row vectors
			% for easier for-loop iteration
			values = struct();
			for param = parameters
				values.(param{:}) = reshape(unique(self.sampleMap.(param{:}), 'stable'), 1, []);
			end
		end
		
		
		function sampleIDs = getSampleIDs(self, treatments)
			% Returns an array of sample IDs corresponding with the given treatments in the order requested
			%
			%	sampleIDs = self.getSampleIDs(treatments)
			%
			%	Inputs
			%		treatments	<struct> A struct where the fields correspond with 
			%					table headers in sampleMap and the values are arrays 
			%					of treatment parameters to be matched in the table.
			%
			%	Outputs
			%		sampleIDs	<numeric> A matrix of sample IDs. Each dimension
			%					in the matrix corresponds with one field of the
			%					'treatments' input. The order of the dimensions
			%					corresponds with the order fields were added to
			%					the struct (since that is the order they pop out
			%					when using the fieldnames() function). 

			% Check treatment requests
			validateattributes(treatments, {'struct'}, {}, mfilename, 'treatments', 1);
			treatmentFields = reshape(fieldnames(treatments), 1, []);
			numTreatments = numel(treatmentFields);
			validTreatments = self.sampleMap.Properties.VariableNames;
			badFields = setdiff(treatmentFields, validTreatments);
			assert(isempty(badFields), ...
				'Field not in sampleMap: %s\n', badFields{:});
			
			% This is used to help rotate the sample IDs into the requested order
			numParams = zeros(1, numTreatments);
			
			% Iterate over fields and order the sample IDs based on matching
			% treatment parameters in the order treatments are requested.
			for i = 1:numTreatments
				
				f = treatmentFields{i};
				treatmentParams = reshape(treatments.(f), 1, []); % Turn into row vector so we can do parallel comparisons
				if ischar(treatmentParams), treatmentParams = {treatmentParams}; end % Force cell array for consistency
				numParams(i) = numel(treatmentParams);
				
				% Extract logical indexes
				matchingSamples = ismember(self.sampleMap.(f), treatmentParams); 
				
				if (i == 1)
					% On first treatment, extract IDs exactly
					IDs = matchingSamples;
				else
					% For subsequent treatments, we add a dimension to the index
					% array which we can use to rapidly and easily find the
					% requested samples in treatment order.
					
					% The previous IDs are first replicated into the next
					% highest dimension (N)
% 					N = ndims(IDs) + 1;
% 					repDims = ones(1, N);
% 					repDims(end) = numParams(i);
% 					IDs = repmat(IDs, repDims);
% 					
% 					% Next we take the current matchingSamples and permute them
% 					% so that their columns now extend into the new dimension
% 					matchingSamples = permute(matchingSamples, [1, N : -1 : 2]);
% 					
% 					% Now we have to replicate matchingSamples into all
% 					% dimensions between 1 and N
% 					sizeIDs = size(IDs); % Form: [numSamples, [prevParams], currParam]
% 					matchingSamples = repmat(matchingSamples, [1, sizeIDs(2:end-1), 1]);
% 					
					IDs = (IDs & matchingSamples);
				end
			end
% 			whos IDs
			% Extract sample IDs as numbers
			[linearSampleIDs, ~] = find(IDs);
% 			numel(linearSampleIDs)
% 			max(linearSampleIDs)
			if (numel(numParams) > 1)
				sampleIDs = reshape(linearSampleIDs, numParams);
			else
				sampleIDs = linearSampleIDs;
			end
		end
		
		
		function dataMatrix = slice(self, sampleIDs, sliceParams)
			% Slices the data, returning an N x M matrix of data from N cells in M channels. 
			%
			%	dataMatrix = self.slice(sampleID, sliceParams)
			%
			%	Inputs
			%		sampleIDs		<integer> The sample(s) to slice as given by
			%						the numerical sample ID(s).
			%		sliceParams		<struct> Optional, struct with optional fields:
			%						'channels': <cell, char>, defaults to self.channels
			%						'dataType': <char>, defaults to 'raw'
			%						'gate':		<char>, defaults to no gate
			%						'bins':		<numeric> defaults to all cells
			%									An Nx1 set of numerical bin IDs or an 
			%									NxD set of bin coordinates where D = #
			%									of bin channels. 
			%									Automatically forces 'gate' and 'dataType' 
			%									to be 'self.binGate' and 'self.binDataType'
			%									regardless of whether they are given or not
			%
			%	Ouputs
			%		dataMatrix		<double> N x M matrix of data from the given
			%						sample where N is the number of cells in the
			%						returned data and M is the number of
			%						channels requested. 
			
			% Check and extract inputs
			[sliceChannels, sliceDataType, sliceGates] = checkInputs_slice(self);
			
			% Slice out data
			dataMatrix = [];
			for s = 1:numel(sampleIDs)
				sID = sampleIDs(s);
				dataS = zeros(sum(sliceGates{s}), numel(sliceChannels));
				for ch = 1:numel(sliceChannels)
					dataS(:, ch) = self.sampleData(sID).(sliceChannels{ch}).(sliceDataType)(sliceGates{s});
				end
				dataMatrix = [dataMatrix; dataS];
			end
			
			
			 % --- Helper Functions --- %
			
			
			function [sliceChannels, sliceDataType, sliceGates] = checkInputs_slice(self)
				% Validates slice properteis and that the sampleID is valid
				
				% Ensure sampleID is a valid integer
				validateattributes(sampleIDs, {'numeric'}, {'positive'}, ...
					mfilename, 'sampleIDs', 1);
				sampleIDs = unique(round(sampleIDs)); 
				assert(all(sampleIDs <= numel(self.sampleData)), ...
					'At least one sampleID is too large!')
				
				% Check and update slice parameters as needed
				if exist('sliceParams', 'var')
					validateattributes(sliceParams, {'struct'}, {}, mfilename, 'sliceParams', 2);
				else
					sliceParams = struct();
				end
				
				if ismember('channels', fieldnames(sliceParams))
					sliceChannels = sliceParams.channels;
					if ischar(sliceChannels), sliceChannels = {sliceChannels}; end % Force cell
					sliceChannels = reshape(sliceChannels, 1, []); % Force row vector
					badChannels = setdiff(sliceChannels, self.channels);
					assert(isempty(badChannels), ...
							'Channel not allowed: %s\n', badChannels{:});
				else
					sliceChannels = self.channels; % Default is all channels
				end
				
				if ismember('dataType', fieldnames(sliceParams))
					validatestring(sliceParams.dataType, self.dataTypes, mfilename, 'sliceParams.dataType');
					sliceDataType = sliceParams.dataType;
				else
					sliceDataType = 'raw'; % Default is raw data
				end
				
				if ismember('gate', fieldnames(sliceParams))
					validatestring(sliceParams.gate, self.gateNames, mfilename, 'sliceParams.gate');
					sliceGates = cell(1, numel(sampleIDs));
					for id = 1:numel(sampleIDs)
						sampleID = sampleIDs(id);
						sliceGates{id} = self.sampleData(sampleID).gates.(sliceParams.gate);
					end
				else
					for id = 1:numel(sampleIDs)
						sampleID = sampleIDs(id);
						sliceGates{id} = true(self.numCells(sampleID), 1); % Default is all cells
					end
				end
				
				if (ismember('bins', fieldnames(sliceParams)) && ~isempty(self.bins))
					validateattributes(sliceParams.bins, {'numeric'}, {}, mfilename, 'sliceParams.bins');
					assert(size(sliceParams.bins, 2) == 1 || ...
						   size(sliceParams.bins, 2) == numel(fieldnames(self.binInputs)), ...
						   'Bin IDs formatted incorrectly!\n')
					assert(size(sliceParams.bins, 1) <= numel(self.bins{sampleID}), ...
						   'Too many bins requested!\n')
					
					% Convert bin indexes to linear 
					if (size(sliceParams.bins, 2) == 1)
						binIdxs = sliceParams.bins';
					else
						spBins = num2cell(sliceParams.bins); % Must convert to cell for dynamic function input
						binIdxs = zeros(size(spBins, 1), 1);
						for bi = 1:numel(binIdxs)
							binIdxs(bi) = sub2ind(size(self.bins{sampleIDs(1)}), spBins{bi, :});
						end
					end
					
					% Override sliceDataType with bin varieties so that the bin 
					% subsampling works right. No need to with gate since the
					% binning doesn't *depend* on the gate. 
					sliceDataType = self.binDataType;
					
					% Extract cells in the requested bins
					for id = 1:numel(sampleIDs)
						sampleID = sampleIDs(id);
						cellsInBins = [];
						for b = binIdxs % Sequentially looks at each bin ID
							cellsInBins = [cellsInBins, self.bins{sampleID}{b}];
						end
						
						inBin = false(size(sliceGates{id}));  % Make 'gate' for cells in bins
						inBin(cellsInBins) = true;		 % Fill out logical index array
						sliceGates{id} = (sliceGates{id} & inBin); % Combine w/ sliceGate for simplicity
					end
				end
			end
		end
		
		
		function threshGate(self, channels, mode, thresh)
			% Thresholds cells in the given channel(s) using the given mode to
			% select thresholded cells from multiple channels. 
			%
			% Individual gates are generated for each given channel, and are
			% then combined together with FlowData.crossGates() using the
			% provided mode. 
			%
			% Uses raw data for thresholding, no matter what stage of processing
			% the data is in. Since the gate is just a logical index array, it
			% generally doesn't matter which dataType is used (though
			% compensation does have a slight effect on distributions).
			%
			%	self.threshGate(channels, mode)
			%
			%	Inputs
			%		channels	<char, cell> A cell list of channels to threshold on. 
			%					(A single string for one channel is also accepted)
			%		mode		<char> {'or', 'and'} - determines how thresh gates are crossed
			%		threshVal	(optional) Either a specific threshold value to use or 
			%					the keyword 'auto', which automatically determines
			%					thresholds from the wtData 99.9th %ile
			%					If numeric, the thresh value given should be a
			%					true (non-logical, non-MEF/MEFL) number. A
			%					single number can be given to threshold in each
			%					channel, or individual numbers for each channel 
			%					can be given. 
			
			checkInputs_threshGate(self);
			
			thrGateNames = cell(1, numel(channels));
			for ch = 1:numel(channels)
				chan = channels{ch};
				
				if (ischar(thresh) && strcmpi(thresh, 'auto'))
					% Thresh value is found using the untransfected (wt) cells
					% --> Set at 99.9th %ile of fluorescence in a channel
					threshVal = prctile(self.controlData(end).(chan).raw, 99.9);
				elseif isnumeric(thresh)
					% Thresh value given as an option
					threshVal = thresh(ch); % Ignore extra entries
				else
					% Thresh value is found by marking a line on a graph
					combData = [];
					for i = 1:self.numSamples
						combData = [combData; self.sampleData(i).(chan).raw(1:self.numSamples:end)];
					end
					
					figure();
					ax = gca(); hold(ax, 'on')
					histogram(ax, Transforms.lin2logicle(combData))
					title('Draw a line to set an x-axis threshold', 'fontsize', 16)
					ylabel('Count', 'fontsize', 14)
					xlabel(strrep(chan, '_', '-'), 'fontsize', 14)
					Plotting.biexpAxes(ax, true, false);
					
					h = imline();
					position = wait(h);
					
					% Fit the points to get a line, then use that to find the
					% x-intercept (threshold value)
					fit = polyfit(position(:, 1), position(:, 2), 1);
% 					plot(ax, (0:0.5:4.5), fit(1) * (0:0.5:4.5) + fit(2), 'r-');
					threshVal = Transforms.logicle2lin(-fit(2) / fit(1));
				end
				
				% Find cells passing threshold and record gate
				thrGateNames{ch} = ['TH_', self.SHORT_COLORS.(chan)];
				for i = 1:self.numSamples
					passThresh = (self.sampleData(i).(chan).raw >= threshVal);
					self.sampleData(i).gates.(thrGateNames{ch}) = passThresh;
				end
				
				% Do the same for controls
				for i = 1:numel(self.controlData)
					if ~isempty(self.controlData(i).(chan)) % Some tcData will be empty 
						passThresh = (self.controlData(i).(chan).raw >= threshVal);
						self.controlData(i).gates.(thrGateNames{ch}) = passThresh;
					end
				end
			end
			
			self.addGates(thrGateNames);
			if (numel(thrGateNames) > 1) % Not necessary for < 2 channels
				self.crossGates(thrGateNames, mode);
			end
			
			fprintf(1, 'Finished thresholding\n');
			
			
			% --- Helper Functions --- %
			
			
			function checkInputs_threshGate(self)
				validateattributes(channels, {'char', 'cell'}, {}, mfilename, 'gates', 1);
				if ischar(channels), channels = {channels}; end % For simplicity
				badChannels = setdiff(channels, self.channels);
				assert(isempty(badChannels), 'Channel not valid: %s\n', badChannels{:});

				if (numel(channels) > 1) % Not necessary for < 2 channels
					validatestring(mode, {'and', 'or'}, mfilename, 'mode', 2);
				end
				
				if exist('thresh', 'var')
					validateattributes(thresh, {'char', 'numeric'}, {}, mfilename, 'threshVal', 3);
					if isnumeric(thresh)
						if (numel(thresh) == 1)
							thresh = repmat(thresh, size(channels));
						end
						assert(numel(thresh) == numel(channels), ...
							'# threshold values given different than # of channels given');
					end
				else
					thresh = false;
				end
			end
		end
				
		
		function crossGates(self, gates, mode)
			% Crosses the given gates using the given crossing mode
			% A new gate is created with name in the following form: 
			%		'gate1_gate2_gate3_[...]_gateN'
			%
			%	self.crossGates(gates, mode)
			%
			%	Inputs
			%		gates	<cell> A cell list of gate names to cross. Must be >= 2 gates
			%		mode	<char> {'or', 'and'} - determines how gates are crossed
			
			% Check inputs
			validateattributes(gates, {'cell'}, {'vector'}, mfilename, 'gates', 1);
			validatestring(mode, {'and', 'or'}, mfilename, 'mode', 2);
			
			badGates = setdiff(gates, self.gateNames);
			assert(isempty(badGates), 'Gate does not exist: %s\n', badGates{:});
			newGateName = [gates{1}, sprintf('_%s', gates{2:end})];
			
			% Do gate crossing
			switch mode
				case 'or'
					for i = 1:self.numSamples
						self.sampleData(i) = crossOR(self.sampleData(i), newGateName);
					end
					for i = 1:numel(self.controlData)
						if ~isempty(self.controlData(i).(self.channels{1})) % Some tcData will be empty 
							self.controlData(i) = crossOR(self.controlData(i), newGateName);
						end
					end
				case 'and'
					for i = 1:self.numSamples
						self.sampleData(i) = crossAND(self.sampleData(i), newGateName);
					end
					for i = 1:numel(self.controlData)
						if ~isempty(self.controlData(i).(self.channels{1})) % Some tcData will be empty 
							self.controlData(i) = crossAND(self.controlData(i), newGateName);
						end
					end
			end
			
			self.addGates(newGateName);
			fprintf(1, 'Finished crossing gates\n')
			
			
			% --- Helper functions --- %
			
			
			function data = crossOR(data, newGateName)
				% Crosses w/ OR logic
				
				crossedGates = false(size(data.gates.(gates{1})));
				for g = 1:numel(gates)
					crossedGates = (crossedGates | data.gates.(gates{g}));
				end
				data.gates.(newGateName) = crossedGates;
			end
			
			
			function data = crossAND(data, newGateName)
				% Crosses w/ AND logic 
				
				crossedGates = true(size(data.gates.(gates{1})));
				for g = 1:numel(gates)
					crossedGates = (crossedGates & data.gates.(gates{g}));
				end
				data.gates.(newGateName) = crossedGates;
			end
		end
		
		
		function [gatePcts] = printGatePcts(self, options)
			% Prints the percent of cells passing each defined gate
			%
			%	Inputs
			%		options		<cell, char> Optional input arguments
			%						'noPrint'	Do not print output
			%
			%	Outputs
			%		gatePcts	<table> Gate percentages in a table format.
			%					Table headers are gate names. 
			
			% Check inputs
			if exist('options', 'var')
				if ischar(options), options = {options}; end
			else
				options = {};
			end
			
			numGates = numel(self.gateNames);
			gatePcts = zeros(self.numSamples, numGates);
			
			for g = 1:numGates
				gate = self.gateNames{g};
				
				for si = 1:self.numSamples
					gatePcts(si, g) = mean(self.sampleData(si).gates.(gate)) * 100;
				end
			end
			
			% Convert to readable table
			gatePcts = array2table(gatePcts);
			gatePcts.Properties.VariableNames = self.gateNames;
			
			% Append to Sample Map
			gatePcts = [self.sampleMap, gatePcts];
			
			if ~ismember('noPrint', options)
				disp(gatePcts)
			end
		end
	end
	
	
	methods (Access = private)
		
% 		function handlePropEvents(src, event)
% 			% Function for handling changes to observable properties
% 			
% 			self = event.AffectedObject;
% 			
% 			fprintf(1, 'Property changed: %s\n', src.Name);
% 			switch src.Name
% 				case {'binInputs', 'binDataType', 'binGate'}
% 					self.bin(self.binInputs, self.binDataType, self.binGate)
% 				
% 				case {'test'}
% 					fprintf(1, 'Test set to: %s\n', event.AffectedObject.test);
% 			end
% 		end
		
		
		function addDataTypes(self, dataTypes)
			% Adds the given dataTypes to the dataTypes property if they are new
			
			validateattributes(dataTypes, {'cell', 'char'}, {}, mfilename, 'dataTypes', 1);
			
			if ~iscell(dataTypes), dataTypes = {dataTypes}; end % For simplicity
			
			% Add dataTypes
			added = {};
			for dt = 1:numel(dataTypes)
				if ~any(strcmp(dataTypes{dt}, self.dataTypes))
					self.dataTypes = [self.dataTypes, dataTypes(dt)];
					added = [added, dataTypes(dt)];
				end
			end
			
			if ~isempty(added), fprintf(1, 'Added dataType: %s\n', added{:}); end
		end
		
		
		function addGates(self, gates)
			% Adds the given gates to the gateNames property if they are new
			
			validateattributes(gates, {'cell', 'char'}, {}, mfilename, 'gates', 1);
			
			if ~iscell(gates), gates = {gates}; end % For simplicity
			
			% Add dataTypes
			added = {};
			for g = 1:numel(gates)
				if ~any(strcmp(gates{g}, self.gateNames))
					self.gateNames = [self.gateNames, gates(g)];
					added = [added, gates(g)];
				end
			end
			
			if ~isempty(added), fprintf(1, 'Added gate: %s\n', added{:}); end
		end
	end
end