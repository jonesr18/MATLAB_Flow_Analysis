classdef FlowData < handle
	% A data structure for managing flow cytometry data. Uses and references
	% several other classes and functions in the repository. 
	%
	%	Visible Properties
	%
	%		name			<char>		Experiment name
	%		date			<char>		Experiment start date
	%		folder			<char>		The full path name for the experiment
	%		cytometer		<char>		The cytometer used for the experiment
	%
	%		numSamples		 <numeric>	The number of data samples
	%		numCells		 <array>	The number of cells in each sample
	%		sampleData		 <struct>	Sample fluorescence and gate data in standard struct
	%		sampleMap		 <table>	Experimental information for samples
	%		dataTypes		 <cell>		Cell array of data types 
	%									('raw', 'comp', 'mefl', etc)
	%		gateNames		 <cell>		Cell array of gate names (strings)
	%		gatePolygons	 <struct>	Mapping between gate names and polygons for sampleData
	%		channels		 <cell>		Cell array of channel names
	%		controlData		 <struct>	Similar to sampleData but for controls.
	%									The order of controls should be the same
	%									as their corresponding channels. 
	%		controlFolder	 <char>		The full path name for the controls
	%		coefficients	 <numeric>	Compensation coefficients
	%		autfluor		 <numeric>  Autofluorescence per channel
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
	%		numBins			 <numeric>	Number of bins
	%		binSizes		 <numeric>	Number of bins in each dimension (order depends 
	%									on the order of channels in binInputs).
	%		binInputs		 <struct>	Struct w/ bin channels as fields and edges as values
	%		binDataType		 <char>		Binned dataType ('raw', 'mComp', 'mefl', etc)
	%		unnamedOps		 <numeric>	% # Of operate() calls with no newDataType input
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
	%
	% Written By
	% Ross Jones
	% jonesr18@mit.edu
	% Weiss Lab, MIT
	
	%#ok<*AGROW>	
	
	properties (SetObservable, Hidden)
		test = '';
	end
	
	
	properties (SetAccess = private)
		name = '';					% Experiment name
		date = '';					% Experiment start date
		folder = '';				% Experiment full path name
		cytometer = '';				% Experiment cytometer used
		
		numSamples = 0;				% The number of data samples
		numControls = 0;			% The number of controls
		numCells = [];				% The number of cells in each sample
		numCellsControls = [];		% The number of cells in each control
		sampleData = struct();		% Sample fluorescence and gate data in standard struct
		sampleMap = table();		% Experimental information for samples
		
		gateNames = {};				% Cell array of gate names (strings)
		gatePolygons = struct();	% Mapping between gate names and polygons for sampleData
		dataTypes = {'raw'};		% Cell array of data types ('raw', 'mComp', 'mefl', etc)
		channels = {};				% Cell array of channel names
		
		controlData = struct();		% Similar to sampleData but for controls
		controlFolder = '';			% Controls full path name
		coefficients = [];			% Compensation coefficients
		autofluor = [];				% Autofluorescence per channel
		compDataType = '';			% Data type used for compensation
		
		beadFitsControls = table();	% Table containing MEF unit fits for controls
		beadFitsSamples = table();	% Table containing MEF unit fits for samples
		meflConversions = struct(); % Mapping between channel names to MEF-MEFL conversion factors 
		
		bins = {};					% Cell array where each element corresponds with a data sample
		numBins = 0;				% Number of bins
		binSizes = [];				% Number of bins in each dimension
		binInputs = struct();		% Struct w/ bin channels as fields and edges as values
		binDataType = '';			% Binned dataType ('raw', 'mComp', 'mefl', etc)
		
		logicleParams = struct( ... % Logicle transformation parameters to use
			'T', 2^18, ...
			'M', 4.5, ...
			'r', -150, ...
			'MEF', Transforms.MEF_CONVERSION_FACTOR);
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
		
		unnamedOps = 0;					% # Of operate() calls with no newDataType input
	end
	
	
	properties (Access = private, Constant)		
		SHORT_COLORS = struct( ...
				'BUV_396_A', 'V', ...
				'Pacific_Blue_A', 'B', ...
				'FITC_A', 'Y', ...
				'PE_A', 'O', ...
				'PE_YG_A', 'O', ...			% Koch LSRII-HTS2 Name
				'PE_Texas_Red_A', 'R', ...
				'PE_TxRed_YG_A', 'R', ...	% Koch LSRII-HTS2 Name
				'APC_Cy7_A', 'I');
	end
	
	
	methods (Access = public)
		
		function self = FlowData(dataFnames, channels, expDetails)
			% Initializes the FlowData object 
			% by importing data from the given files, which should correspond
			% with the given sample map and contain data in the given channels.
			%
			%	self = FlowData(dataFnames, channels, expDetails)
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
			%		expDetails		<struct> A struct with the following fields
			%						recording experimental details:
			%							name		|	<char>
			%							date		|	<char>
			%							folder		|	<char> (full-path folder name)
			%							sampleMap*	|	<char>
			%							cytometer	|	<char>
			%
			%						*sampleMap: name of the .txt file containing 
			%						treatments information for each sample. 
			%						**See example file in source folder.
			%
			%	Outputs
			%		self			A handle to the object
			
			% When loading, we need to make a blank version of the object
			if strcmpi(dataFnames, 'load'), return, end
			
			[dataStruct, sampleMapFname] = zCheckInputs(self);
			
			% Extract experiment details
			self.date = expDetails.date;
			self.name = expDetails.name;
			self.folder = expDetails.folder;
			self.cytometer = expDetails.cytometer;
			
			% Defaults
			numReplicates = 1;
			IDs = 1:numel(dataStruct);
			
			% Extract sample map
			if ~isempty(sampleMapFname)
				try % Don't want to throw away loaded data
					sampleMap = readtable(sampleMapFname, 'delimiter', '\t');
					assert(height(sampleMap) == numel(dataStruct), ...
						'Sample map (%d) contains incorrect number of rows (%d)', ...
						height(sampleMap), numel(dataStruct))
					
					% "Replicate" column allows data to be automatically combined
					if ismember('Replicate', sampleMap.Properties.VariableNames)
						% Overwrite numReplicates and IDs based on data condensation 
						numReplicates = max(sampleMap.Replicate);
						IDs = zeros(numReplicates, height(sampleMap) / numReplicates);
						for r = 1:numReplicates
							IDs(r, :) = find(sampleMap.Replicate == r);
						end
						% Wean sampleMap to only be have one replicate
						sampleMap = sampleMap(sampleMap.Replicate == 1, :);
					end
					self.sampleMap = sampleMap;
				catch ME
					warning('Error occured while reading sample map, skipping...')
					disp(ME)
				end
			end
			
			% Extract data from the given channels
			self.numSamples = size(IDs, 2);
			self.numCells = zeros(1, self.numSamples);
			for i = 1:size(IDs, 2)
				
				nObs = 0;
				
				for r = 1:numReplicates
					nObs = nObs + dataStruct(IDs(r, i)).nObs;
				end
				self.numCells(i) = nObs;
				
				% Extract desired color channels
				for ch = channels
					sd = [];
					for r = 1:numReplicates
						sd = [sd; dataStruct(IDs(r, i)).(ch{:}).raw];
					end
					self.sampleData(i).(ch{:}).raw = sd;
				end
				self.sampleData(i).nObs = nObs;
				
				% Extract scatter channels
				for ch = Gating.SCATTER_CHANNELS
					sds = [];
					for r = 1:numReplicates
						sds = [sds; dataStruct(IDs(r, i)).(ch{:}).raw];
					end
					self.sampleDataScatter(i).(ch{:}).raw = sds;
				end
				self.sampleDataScatter(i).nObs = nObs;
			end
			self.channels = channels;
			
			% Add listeners for settable public properties
% 			addlistener(self, 'test', 'PostSet', @self.handlePropEvents);
% 			addlistener(self, {'binInputs', 'binDataType'}, ...
% 				'PostSet', @self.handlePropEvents);
			
			fprintf(1, 'Finished constructing FlowData object\n')
			
			
			% -- Helper functions -- %
			
			
			function [dataStruct, sampleMapFname] = zCheckInputs(self)
				validateattributes(dataFnames, {'cell', 'char'}, {}, mfilename, 'dataFilenames', 1);
				validateattributes(channels, {'cell', 'char'}, {}, mfilename, 'channels', 2);
				validateattributes(expDetails, {'struct'}, {}, mfilename, 'exptDetails', 3);
				
				% Convert channels to cell array if single char value is given
				if ischar(channels), channels = {channels}; end
				channels = reshape(channels, 1, []); % Ensure row vector
				
				% Check required experiment details are present
				requiredFields = {'date', 'name', 'folder', 'cytometer'};
				missingFields = setdiff(requiredFields, fieldnames(expDetails));
				assert(isempty(missingFields), 'Experiment details missing field: %s\n', missingFields{:});
				assert(logical(exist(expDetails.folder, 'file')), 'Experiment folder does not exist!');
				if ~(expDetails.folder(end) == filesep)
					expDetails.folder = [expDetails.folder, filesep];
				end
				
				% Check sampleMap is a real file
				if isfield(expDetails, 'sampleMap')
					sampleMapFname = self.convertToFullFile(expDetails.sampleMap, expDetails.folder);
					if ~logical(exist(sampleMapFname, 'file'))
						warning('Sample map not found: %s, skipping...', sampleMapFname)
					end
				else
					sampleMapFname = [];
				end
				
				% Get data filenames
				dataFnames = self.convertToFullFile(dataFnames, expDetails.folder);
				ds = FlowAnalysis.openFiles(dataFnames{1});
				
				% Check channels are present in dataStruct
				badChannels = setdiff(channels, fieldnames(ds));
				assert(isempty(badChannels), ...
					'Channel not in dataStruct: %s\n', badChannels{:});
				
				% Import all data
				if (numel(dataFnames) > 1)
					dataStruct = [ds, FlowAnalysis.openFiles(dataFnames{2:end})];
				else
					dataStruct = ds;
				end
			end
		end
		
		
		function self = addControls(self, controlFolder, wildTypeFname, singleColorFnames, twoColorFnames)
			% Adds wild-type, single-color, and two-color (optional) data to the dataset
			% so that we can do compensation (single-colors) and MEFL conversion (two-colors).
			%
			% The method generates self.controlData, a struct array where single-color  
			% data from channel X is in position X, two-color data from channel X is
			% in position 2*X (if applicable), and wild-type data is in the last
			% position (regardless of the presence of two-color data). 
			% 
			%	self.addControls(controlFolder, wildTypeFname, singleColorFnames, twoColorFnames)
			%
			%	Inputs
			%		controlFolder		<char> The full-path folder name for controls
			%
			%		wildTypeFnames		<cell, char> Wild-type cell data
			%		
			%		singleColorFnames	<cell, char> Single-color controls data files
			%							** Order of colors should coincide with
			%							the order of FlowData.channels
			%
			%		twoColorFnames		(optional) <cell, char> Two-color controls data files
			%							** As with single colors, the order should match
			%							those in self.channels, but with no yellow/green 
			%							file (since all other colors are converted to MEFL 
			%							units using the FITC channel).
			
			[wildTypeData, singleColorData, twoColorData] = zCheckInputs_addControls(self);
			FITC_IDX = find(strcmpi('FITC_A', self.channels));
			self.controlFolder = controlFolder;
			
			% Extract data
			self.controlData = extractData([self.channels, {'nObs'}], ...
						wildTypeData, singleColorData, twoColorData, FITC_IDX);
			
			% Extract scatter data
			self.controlDataScatter = extractData([Gating.SCATTER_CHANNELS, {'nObs'}], ...
						wildTypeData, singleColorData, twoColorData, FITC_IDX);
			
			% Find number of cells per control
			self.numControls = numel(self.controlData);
			self.numCellsControls = zeros(self.numControls, 1);
			for ci = 1:self.numControls
				nc = self.controlData.nObs;
				if isempty(nc), nc = 0; end
				self.numCellsControls = nc;
			end
			
			self.controlsAdded = true;
			fprintf(1, 'Finished adding controls\n');
			
			
			% --- Helper Functions --- %
			
			
			function [wildTypeData, singleColorData, twoColorData] = zCheckInputs_addControls(self)
				
				validateattributes(controlFolder, {'char'}, {}, mfilename, 'controlFolder', 1);
				assert(logical(exist(controlFolder, 'file')), 'Controls folder does not exist!');
				if ~(controlFolder(end) == filesep), controlFolder = [controlFolder, filesep]; end
				
				validateattributes(wildTypeFname, {'cell', 'char'}, {}, mfilename, 'wildTypeFname', 2);
				validateattributes(singleColorFnames, {'cell', 'char'}, {}, mfilename, 'singleColorFnames', 3);
				
				% Convert to cell arrays if necessary for convenience
				if ischar(wildTypeFname), wildTypeFname = {wildTypeFname}; end
				if ischar(singleColorFnames), singleColorFnames = {singleColorFnames}; end
				
				% Check number of scFiles
				assert(numel(singleColorFnames) == numel(self.channels), ...
					'Incorrect number of single color controls');
				
				% Add full-path to filenames
				wildTypeFname = self.convertToFullFile(wildTypeFname, controlFolder);
				singleColorFnames = self.convertToFullFile(singleColorFnames, controlFolder);
				
				% Open files
				wildTypeData = FlowAnalysis.openFiles(wildTypeFname{:});
				singleColorData = FlowAnalysis.openFiles(singleColorFnames{:});
				
				% Add twoColorData if applicable
				if exist('twoColorFnames', 'var')
					validateattributes(twoColorFnames, {'cell', 'char'}, {}, mfilename, 'twoColorFnames', 4);
					if ischar(twoColorFnames), twoColorFnames = {twoColorFnames}; end
					assert(numel(twoColorFnames) == sum(~strcmpi('FITC_A', self.channels)), ...
						'Incorrect number of two color controls');
					twoColorFnames = self.convertToFullFile(twoColorFnames, controlFolder);
					twoColorData = FlowAnalysis.openFiles(twoColorFnames{:});
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
		
		
		function self = gate(self, onlyP1)
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
			onlyP1 = (exist('onlyP1', 'var') && all(logical(onlyP1)));
			self.onlyP1 = onlyP1;
			if strcmpi(self.cytometer, 'Koch-LSRII-HTS2'), swap = true; else, swap = false; end
			
			% Setup new directories for gates
			gateDirControls = [self.controlFolder, 'Gating', filesep];
			if ~exist(gateDirControls, 'file')
				mkdir(gateDirControls)
			end
			gateDirSamples = [self.folder, 'Gating', filesep];
			if ~exist(gateDirSamples, 'file')
				mkdir(gateDirSamples)
			end
			
			% Check if gates have already been made for this data first, then
			% process if necessary. 
			gatesSaveName = [self.date, '_', self.name];
			
			gatesFnameControls = [gateDirControls, 'Controls_GatePolygons.mat'];
			if exist(gatesFnameControls, 'file')
				% Load existing control gates
				load(gatesFnameControls, 'gateP1c', 'gateP2c', 'gateP3c');
			else
				% Do manual control gating
				[gateP1c, gateP2c, gateP3c, gateFigs] = Gating.standardGating(self.controlDataScatter, onlyP1, swap);
				save(gatesFnameControls, 'gateP1c', 'gateP2c', 'gateP3c');
				for f = fieldnames(gateFigs)'
					saveas(gateFigs.(f{:}), [gateDirControls, 'Controls_gate', f{:}, 'c']);
				end
			end
			
			gatesFnameSamples = [gateDirSamples, gatesSaveName '_GatePolygons.mat'];
			if exist(gatesFnameSamples, 'file')
				% Load existing sample gates
				load(gatesFnameSamples, 'gateP1s', 'gateP2s', 'gateP3s');
			else
				% Do manual sample gating
				[gateP1s, gateP2s, gateP3s, gateFigs] = Gating.standardGating(self.sampleDataScatter, onlyP1, swap);
				save(gatesFnameSamples, 'gateP1s', 'gateP2s', 'gateP3s');
				for f = fieldnames(gateFigs)'
					saveas(gateFigs.(f{:}), [gateDirSamples, gatesSaveName, '_gate', f{:}, 's']);
				end
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
			
			% Transfer gate logicals to externally accesible data
			for cd = 1:self.numControls
				self.controlData(cd).gates = self.controlDataScatter(cd).gates;
			end
			for sd = 1:self.numSamples
				self.sampleData(sd).gates = self.sampleDataScatter(sd).gates;
			end
			
			fprintf(1, 'Finished standard gating\n')
		end
		
		
		function self = customGate(self, gateName, sampleIDs, sliceParams, axScale)
			% Creates a custom gate using data from the given samples 
			% and the given slice parameters.
			%
			%	self.customGate(gateName, sampleIDs, sliceParams)
			%
			%	Inputs
			%
			%		gateName		<char> The name of the new gate
			%
			%		sampleIDs		<integer> The sample(s) to slice as given by
			%						the numerical sample ID(s).
			%							(Optional, defaults to all cells)
			%
			%		sliceParams		<struct> Optional, struct with optional fields:
			%						'channels':  <cell, char>, defaults to self.channels
			%									 (Only first 2 channels are used)
			%						'dataType':  <char>, defaults to 'raw'
			%						'gate':		 <char>, defaults to no gate
			%						'equalize':  <logical>, TRUE returns an equal
			%									 number of points from each sample
			%									 (default = FALSE), operation
			%									 performed before applying bins
			%						'numPoints': The number of cells to extract
			%									 per sample (if the minimum number 
			%									 of cells per sample is lower, the 
			%									 method will use that value). 
			%						'controls':  <logical> TRUE slices from
			%									 controlData rather than sampleData
			%									 (default = FALSE)
			%						'bins':		 <numeric> defaults to all cells
			%									 An Nx1 set of numerical bin IDs or an 
			%									 NxD set of bin coordinates where D = #
			%									 of bin channels. 
			%									 Automatically forces 'dataType' to be 
			%									 'self.binDataType' regardless of whether 
			%									 they are given or not
			%									 <Can input 'all' to select all bins>
			%
			%		axScale			<char> Optional: The axis scaling to use:
			%							'loglog', 'semilogy', 'semilogx',
			%							'linear' (default)
			
			% Check inputs
			zCheckInputs_customGate()
			
			gateDirSamples = [self.folder, 'Gating', filesep];
			if ~exist(gateDirSamples, 'file')
				mkdir(gateDirSamples)
			end
			
			% Extract data and do gating
			outData = self.slice(sampleIDs, sliceParams);
			[~, gatePolygon, gateFig] = Gating.gatePolygon( ...
					outData(:, 1), outData(:, 2), axScale);
			
			% Store gate info
			if (isfield(sliceParams, 'controls') && sliceParams.controls)
				for ci = 1:self.numControls
					if isempty(self.controlData(ci).(self.channels{1}))
						continue % Handle empty tcData
					end
					% Extract sample data and gate
					[outData, sliceGates] = self.slice(ci, sliceParams);
					inGate = Gating.gatePolygon(outData(:, 1), outData(:, 2), axScale, gatePolygon);
					
					% Adjust index based on input gate
					inGateFixed = Gating.fixGateIdxs(sliceGates{1}, inGate);
					self.controlData(ci).gates.(gateName) = inGateFixed;
				end
			else
				for si = 1:numel(self.sampleData)
					% Extract sample data and gate
					[outData, sliceGates] = self.slice(si, sliceParams);
					inGate = Gating.gatePolygon(outData(:, 1), outData(:, 2), axScale, gatePolygon);
					
					% Adjust index based on input gate
					inGateFixed = Gating.fixGateIdxs(sliceGates{1}, inGate);
					self.sampleData(si).gates.(gateName) = inGateFixed;
				end
			end
			
			% Todo adjust for existing gate idxs
			
			% Save gate polygon
			self.gatePolygons.(gateName) = gatePolygon;
			self.addGates(gateName);
				
			% Save gate figure
			if ~isempty(gateFig)
				gateFigFname = [gateDirSamples, self.date, '_', self.name, '_gate', gateName, '.fig'];
				saveas(gateFig, gateFigFname);
			end
			
			
			% --- Helper Functions --- % 
			
			
			function zCheckInputs_customGate()
				
				if exist('sliceParams', 'var')
					validateattributes(sliceParams, {'struct'}, {}, mfilename, 'sliceParams', 2);
				else
					sliceParams = struct();
				end
				
				if ~exist('axScale', 'var')
					axScale = 'linear';
				end
			end
		end
		
		
		function self = convertToMEF(self, beadsControls, beadsSamples, options)
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
			%				   ({'nneg', 'P1_nneg'} if self.onlyP1 = TRUE)
			%
			%	Inputs
			%
			%		beadsControls	<struct> A struct with the following fields:
			%			filename		The .fcs file containing bead data
			%							corresponding with this experiment
			%			type			The name of the type of bead (eg 'RCP-30-5A')
			%			lot				The bead production lot (eg 'AH01')
			%
			%		beadsSamples	<struct> A struct with the same fields as above,
			%						but for the samples rather than controls
			%
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
			zCheckInputs_convertToMEF(self);
			
			% Setup new directory for fitting files/figs
			beadDirControls = [self.controlFolder, 'Calibration', filesep];
			if ~exist(beadDirControls, 'file')
				mkdir(beadDirControls)
			end
			beadDirSamples = [self.folder, 'Calibration', filesep];
			if ~exist(beadDirSamples, 'file')
				mkdir(beadDirSamples)
			end
			
			% Extract MEF units from Transforms class for naming
			MEF_units = Transforms.getBeadUnits(self.channels);
			
			beads = {beadsControls, beadsSamples};
			beadDirs = {beadDirControls, beadDirSamples};
			for b = 1:numel(beads)
				
				% This is unnessesary - remove after ensuring so
% 				if isequaln(beadsControls, beadsSamples)
% 					% If the bead properties are the same, then skip the fitting for
% 					% samples' beads and just use the controls' beads. This is for 
% 					% the case where samples and controls are run the same day.
% 					fitsSamples = fitsControls;
% 					fprintf(1, 'Sample beads match control beads...using them for calibration\n');
% 					break
% 				end
				
				% Setup filenames, check if already exists
				beadSaveName = [beads{b}.date, '_', beads{b}.type, '_', ...
						beads{b}.lot, '_', beads{b}.cytometer];
				mefFname = [beadDirs{b}, beadSaveName '_MEF_Fits.mat'];
				
				if exist(mefFname, 'file')
					fprintf(1, 'Loading pre-computed bead fits\n');
					load(mefFname, 'mefFits');
				else
					% Get MEF fits 
					[mefFits, figFits] = Transforms.calibrateMEF(beads{b}, self.channels, options);

					% Save channel fits (and figures if applicable)
					save(mefFname, 'mefFits')

					% Save figures to bead directory
					if ~isempty(fieldnames(figFits))
						for chID = 1:numel(self.channels)
							saveas(figFits.(self.channels{chID}), ...
									[beadDirs{b}, beadSaveName, '_', self.channels{chID}, '_', MEF_units{chID}, '_Fit.fig']); 
						end
						if isfield(figFits, 'manualPeaks')
							saveas(figFits.manualPeaks, [beadDirs{b}, beadSaveName, '_manualPeaks.fig']);
						end
					end
				end
				
				% Extract fits to unique local variables
				if (b == 1)
					fitsControls = mefFits; 
				else
					fitsSamples = mefFits; 
				end
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
			
			
			function zCheckInputs_convertToMEF(self)
				
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
		
		
		function self = convertToMEFL(self, showPlots)
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
						
			% Check if conversions already exist
			beadDir = [self.controlFolder, 'Calibration', filesep];
			if ~exist(beadDir, 'file')
				error('Controls bead directory not found! It should be set up during MEF calibration')
			end
			meflFname = [beadDir, 'MEFL_Conversions.mat'];
			
			if exist(meflFname, 'file')
				fprintf(1, 'Loading pre-computed MEFL conversions\n');
				load(meflFname, 'meflFits')
			else
				% Compute mefl conversions
				tcData = self.controlData(numel(self.channels) + 1 : 2 * numel(self.channels));
				[meflFits, figFits] = Transforms.calibrateMEFL(tcData, self.channels, 'mef', showPlots);
			
				save(meflFname, 'meflFits');
				
				if ~isempty(fieldnames(figFits))
					for chID = 1:numel(self.channels)
						if ~strcmpi(self.channels{chID}, 'FITC_A')
							saveas(figFits.(self.channels{chID}), ...
									[beadDir, self.channels{chID}, '_MEFL_Conversion.fig']); 
						end
					end
				end
			end
			
			% Add converted data as 'mefl' data type
			for ch = self.channels
				% Add MEFLs for controls
				for i = 1:self.numControls
					if isempty(self.controlData(i).(ch{:})), continue, end % Some tcData will be empty 
					self.controlData(i).(ch{:}).mefl = self.controlData(i).(ch{:}).mef * meflFits.(ch{:});
				end
				
				% Add MEFLs for sample
				for i = 1:self.numSamples
					self.sampleData(i).(ch{:}).mefl = self.sampleData(i).(ch{:}).mef * meflFits.(ch{:});
				end
			end
			
			self.meflConversions = meflFits; 
			self.addDataTypes('mefl');
			self.meflConverted = true;
			fprintf(1, 'Finished converting to MEFL\n');
		end
		
		
		function self = compensate(self, dataType, gates, options)
			% Applies autofluorescence subtraction and matrix-based compensation
			%
			%	self.compensate(dataType, gates, options)
			%
			%	This function adds new dataTypes: {'comp', 'afs'}
			%
			%	Inputs
			%
			%		dataType	<char> Indicates which data type to use.
			%					Can be any dataType in self.dataTypes
			%
			%		gates		<char, cell> Indicates which gate(s) to use.
			%					Can be any gate in self.gateNames, but
			%					preferrably 'P1', or 'P3' if onlyP1 = false
			%					 -> Also accepts a cell array of gate names,
			%						each applied to a control (same order)
			%
			%		options		<struct> (optional) Optional property-value pairs:
			%			'plotsOn':		If TRUE, shows the compensation plots, which
			%							are then saved in the Compensation folder.
			%			'plotLin':		If TRUE (and plotsOn = TRUE), then the
			%							fits are plotted in linear space, rather
			%							than biexponential
			%			'minFunc':		A user-defined function for residual 
			%							minimzation during fitting. 
			%								Default = @(x) sum(x.^2)
			%								(least-squares approximation)
			%			'recompute':	If TRUE, forces re-calculation of the 
			%							compensation fits
			%			'save':			If TRUE, saves the compensation fits and
			%							(if applicable) the generated figures
			%								Default = true
			
			% Check inputs
			zCheckInputs_compensate(self);
			
			% Setup new directory for compensation fits/figs
			compDir = [self.controlFolder, 'Compensation', filesep];
			if ~exist(compDir, 'file')
				mkdir(compDir)
			end
			
			fitFigs = struct();
			
			% Check if compensation has already been done for this data
			% --> If so, load and apply to data!
			compFname = [compDir, 'CompFits.mat'];
			if (exist(compFname, 'file') && ~all(logical(options.recompute)))
				% Load existing sample gates
				fprintf(1, 'Loading pre-computed coefficients and autofluorescence!\n');
				load(compFname, 'coeffs', 'ints');
			else
				% Compute slopes/intercepts for each pair of scData
				scData = cell(1, numel(self.channels));
				for sc = 1:numel(scData)
					
					% This slice works because controlData 1-C are the
					% single-color controls for channels 1-C
					scData{sc} = self.slice(sc, struct( ...
						'controls', true, ...
						'dataType', dataType, ...
						'gate', gates{sc}));
					
					% Code for binnind data first
% 					maxVal = max(slicedData(:));
% 					edges = cell(1, numel(scData));
% 					edges{sc} = linspace(0, maxVal, 50);
% 					edges(setdiff(1:numel(scData), sc)) = {linspace(0, maxVal, 2)};
% 					binIdxs = reshape(FlowAnalysis.simpleBin(slicedData, edges), 1, []);
% 					
% 					binMedians = cell2mat( ...
% 						cellfun(@(x) median(slicedData(x, :), 1), binIdxs, ...
% 						'uniformoutput', false)');
% 					scData{sc} = binMedians(~isnan(binMedians(:, 1)), :); % Remove NaNs
				end
					[coeffs, ints, fitFigs.pre] = Compensation.computeCoeffs( ...
						scData, self.channels, options); 
			end
			
			% Subtract autofluorescence
			self.operate('add', -ints, dataType, 'afs');
			
			% Compensate
			for cd = 1:self.numControls
				if isempty(self.controlData(cd).(self.channels{1})), continue, end
				
				dataMatrix = self.slice(cd, struct( ...
					'controls', true, ...
					'dataType', 'afs'));
				fixedMatrix = Compensation.matrixComp(dataMatrix', coeffs)';
				for ch = 1:numel(self.channels)
					self.controlData(cd).(self.channels{ch}).comp = fixedMatrix(:, ch);
				end
			end
			for sd = 1:self.numSamples
				dataMatrix = self.slice(sd, struct( ...
					'dataType', 'afs'));
				fixedMatrix = Compensation.matrixComp(dataMatrix', coeffs)';
				for ch = 1:numel(self.channels)
					self.sampleData(sd).(self.channels{ch}).comp = fixedMatrix(:, ch);
				end
			end
			
			self.addDataTypes('comp');
			
			% Check compensation
			if options.plotsOn
				scData = cell(1, numel(self.channels));
				for sc = 1:numel(scData)
					% This slice works because controlData 1-C are the
					% single-color controls for channels 1-C
					scData{sc} = self.slice(sc, struct( ...
						'controls', true, ...
						'dataType', 'comp', ...
						'gate', gates{sc}));
					
					% Code for binning data first
% 					maxVal = max(slicedData(:));
% 					edges = cell(1, numel(scData));
% 					edges{sc} = linspace(0, maxVal, 50);
% 					edges(setdiff(1:numel(scData), sc)) = {linspace(0, maxVal, 2)};
% 					binIdxs = reshape(FlowAnalysis.simpleBin(slicedData, edges), 1, []);
% 					
% 					binMedians = cell2mat( ...
% 						cellfun(@(x) median(slicedData(x, :), 1), binIdxs, ...
% 						'uniformoutput', false)');
% 					scData{sc} = binMedians(~isnan(binMedians(:, 1)), :); % Remove NaNs
				end
				[~, ~, fitFigs.post] = Compensation.computeCoeffs( ...
						scData, self.channels, options);
			end
			
			% Save data/figures (default = do saving)
			if (~isfield(options, 'save') || options.save)
				save(compFname, 'coeffs', 'ints')
				if (isfield(fitFigs, 'pre') && ~isempty(fitFigs.pre)) 
					% Doesn't run if data was re-loaded or plots were not requested
					saveas(fitFigs.pre, [compDir, 'pre-comp.fig'])
					saveas(fitFigs.post, [compDir, 'post-comp.fig'])
				end
			end
			
			self.coefficients = coeffs; 
			self.autofluor = ints;
% 			self.compensated = true;
			self.compDataType = dataType;
			fprintf(1, 'Finished compensation\n')
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_compensate(self)
				validatestring(dataType, self.dataTypes, mfilename, 'dataType', 1);
				if ischar(gates), gates = {gates}; end % For simplicity
				for gi = 1:numel(gates)
					validatestring(gates{gi}, self.gateNames, mfilename, 'gate', 2);
				end
				% If one gate given, extend gates to match # of channels
				if (numel(gates) == 1), gates = repmat(gates, 1, numel(self.channels)); end
				
				if ~exist('options', 'var'), options = struct(); end
				if ~isfield(options, 'plotsOn'), options.plotsOn = false; end
				options.plotsOn = all(logical(options.plotsOn));
				if isfield(options, 'minFunc')
					validateattributes(options.minFunc, {'function_handle'}, {}, mfilename, 'options.minFunc', 3)
				else
					options.minFunc = @(x) sum(x.^2);
				end
				if ~isfield(options, 'recompute'), options.recompute = false; end
				options.logicle = self.logicleParams;
				options.doMEF = ismember(dataType, {'mef', 'mefl'}); 
			end
			
		end
		
		
		function self = bin(self, binInputs, binDataType, doPar)
			% Sorts the sample data using the given set of channels/edges into a
			% number of bins using the given dataType for assignments. 
			%
			%	self.bin(binInputs, binDataType)
			%
			%	Inputs
			%		binInputs		<struct> A struct with channel names as keys and 
			%						bin edges as values. The channel names must match 
			%						a field in the data struct with the subfield 'raw'. 
			%						The struct tells the function which channels to bin 
			%						on and where to draw the bins in each dimension. 
			%		
			%		binDataType		<char> The cell dataType to use (eg 'mefl', 'comp')
			%
			%		doPar			<logical> (Optional) Flag to use parallel
			%						computing for binning. Default = false. Runs
			%						with default MATLAB pool generation. 
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
			%		The reason to do this relatively slow binning method is so
			%		that accessing which cells are in each bin operates in O(1)
			%		time, as we can directly query a bin to see which cells are
			%		inside of it, making all operations _after_ binning much 
			%		faster, overall improving the speed of analysis. 
			
			[binChannels, binEdges] = zCheckInputs_bin(self);
			
			slicedData = cell(1, self.numSamples);
			for i = 1:self.numSamples
								
				slicedData{i} = self.slice(i, struct( ...
						'channels', {binChannels}, ...
						'dataType', binDataType, ...
						'equalize', false));
				
			end
			
			binnedData = cell(size(slicedData));
			if doPar
				parfor i = 1:numel(slicedData)
					binnedData{i} = FlowAnalysis.simpleBin(slicedData{i}, binEdges);
				end
			else
				for i = 1:numel(slicedData)
					binnedData{i} = FlowAnalysis.simpleBin(slicedData{i}, binEdges);
				end
			end
			
			self.bins = binnedData;
			self.numBins = numel(binnedData{1});
			bSizes = size(binnedData{1});
			self.binSizes = bSizes(1:numel(binChannels)); % Handles single-channel binning so only one size is given
			self.binInputs = binInputs;
			self.binDataType = binDataType;
			self.binned = true;
			fprintf(1, 'Finished binning\n')
			
			
			% --- Helper Functions --- %
			
			
			function [binChannels, binEdges] = zCheckInputs_bin(self)
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
				
				doPar = exist('doPar', 'var') && all(logical(doPar));
			end
		end
		
		
		function binStats = computeBinStats(self, sampleIDs, sliceParams, metrics)
			% Computes the given bin statistic metrics on the combined data 
			% from the given sampleIDs.
			%
			%	Inputs
			%		sampleIDs		<integer> The sample(s) as given by the numerical 
			%						sample ID(s) to combine and compute bin stats with.
			%
			%		sliceParams		<struct> Optional, struct with optional fields:
			%						'channels': <cell, char>, defaults to self.channels
			%						'dataType': <char>, defaults to self.binDataType
			%						'gate':		<char>, defaults to no gate
			%						'equalize': <logical>, TRUE returns an equal
			%									number of points from each sample
			%									(default = FALSE), operation
			%									performed before applying bins
			%						'numPoints': The number of cells to extract
			%									 per sample (if the minimum number 
			%									 of cells per sample is lower, the 
			%									 method will use that value). 
			%						'bins':		<numeric> defaults to all bins
			%									An Nx1 set of numerical bin IDs or an 
			%									NxD set of bin coordinates where D = #
			%									of bin channels. 
			%									Automatically forces 'dataType' to be 
			%									'self.binDataType' regardless of whether 
			%									they are given or not
			%									Can input 'all' to select all bins.
			%									** Requesting anything but all bins
			%									will cause the output to be 2D rather 
			%									than N+1 D as the data will not be 
			%									reshaped to match the bin sizes.  
			%
			%		metrics			<char, cell> (Optional) A list of metrics to compute
			%						Valid metrics:
			%							'numCells', pctiles: 'p10', 'p50'/'median', 'p90',
			%							'mean', 'geomean', 'stdev', 'geostdev', 'sem', semb*
			%							 * semb = bootstrapped SEM
			%							 -> If no metrics input given, then all
			%								metrics except 'semb' are computed
			%
			%	Outputs
			%		binStats	<struct> A struct where each field is the name of 
			%					a statistical metric and the value is a N+1 D
			%					matrix of the metric values in N bin dimensions 
			%					across C channels (given in sliceParams) in the 
			%					final dimension. 
			%					--> In the case of requesting stats from specific bins
			%						rather than all of them, we do not structure the 
			%						bin stats into their "true" shape and instead return 
			%						a BxC matrix of metric values in each of the 
			%						requested B bins in all the C channels. 
			
			zCheckInputs_computeBinStats(self);
			
			% Set up binStats struct
			sliceBins = sliceParams.bins;
			binAll = (numel(sliceBins) == self.numBins);
			if binAll
				binStatSize = self.binSizes;
			else
				binStatSize = size(sliceBins, 1);
			end
			for m = metrics
				if binAll
					binStats.(m{:}) = zeros([binStatSize, numel(sliceParams.channels)]);
				end
			end
			binStats.numCells = zeros(binStatSize);
			
			% Iterate over bins, compute stats for each one independently
			for b = 1:size(sliceBins, 1)
				
				% Slice out data to compute stats with
				sliceParams.bins = sliceBins(b, :);
				slicedData = self.slice(sampleIDs, sliceParams);
				binStats.numCells(b) = size(slicedData, 1);
				
				for ch = 1:numel(sliceParams.channels)
					
					% Compute linear index
					binStatIdx = b + numel(sliceBins) * (ch - 1);
					
					% Extract channel data
					dataInBin = slicedData(:, ch);
					dataInBin = dataInBin(~isnan(dataInBin));
					posData = (dataInBin > 0);
					
					% Calculate requested metrics
					prctiles = prctile(dataInBin, [10, 50, 90]);
					if ismember('p10', metrics), binStats.p10(binStatIdx) = prctiles(1); end
					if any(ismember({'p50', 'median'}, metrics)), binStats.p50(binStatIdx) = prctiles(2); end
					if ismember('p90', metrics), binStats.p90(binStatIdx) = prctiles(3); end
					if ismember('mean', metrics), binStats.mean(binStatIdx) = mean(dataInBin); end
					if ismember('geomean', metrics), binStats.geomean(binStatIdx) = geomean(dataInBin(posData)); end
					if ismember('stdev', metrics), binStats.stdev(binStatIdx) = std(dataInBin); end
					if ismember('geostdev', metrics), binStats.geostdev(binStatIdx) = geostd(dataInBin(posData)); end
					if ismember('sem', metrics), binStats.sem(binStatIdx) = std(dataInBin) / sqrt(numel(dataInBin)); end
					if ismember('semb', metrics), binStats.semb(binStatIdx) = semBootstrap(dataInBin); end % NOTE: SUPER SLOW
% 					if ismember('ci95', metrics), binStats.CI95(binStatIdx) = ci95(dataInBin); end
				end
			end
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_computeBinStats(self)
				
				validateattributes(sampleIDs, {'numeric'}, {'positive'}, mfilename, 'sampleIDs', 1);
				sampleIDs = unique(ceil(sampleIDs));
				assert(max(sampleIDs) <= self.numSamples, 'Sample ID(s) too large!')
				
				% Only need to check sliceParams for fields that are used by
				% this function, so we let slice() itself check the rest.
				validateattributes(sliceParams, {'struct'}, {}, mfilename, 'sliceParams', 2);
				if isfield(sliceParams, 'channels') 
					if ischar(sliceParams.channels)
						sliceParams.channels = {sliceParams.channels}; % Needed in main function
					end
				else
					sliceParams.channels = self.channels;
				end
				
				% The default behavior is slightly different than slice - we
				% take all bins instead of all cells when no input is given
				if (~isfield(sliceParams, 'bins') || isempty(sliceParams.bins))
					sliceParams.bins = 'all';
				end
				if ischar(sliceParams.bins)
					if strcmpi(sliceParams.bins, 'all')
						sliceParams.bins = (1:self.numBins)';
					else
						error('Bin char input not recognized: %s', sliceParams.bins)
					end
				end
				
				validMetrics = {'numCells', 'p10', 'p50', 'p90', 'median', 'mean', 'geomean', 'stdev', 'geostdev', 'sem', 'semb'};
				if exist('metrics', 'var')
					validateattributes(metrics, {'cell', 'char'}, {}, mfilename, 'metrics', 3);
					if ischar(metrics), metrics = {metrics}; end % For simplicity
					metrics = reshape(metrics, 1, []); % Force row vector
					badMetrics = setdiff(metrics, validMetrics);
					assert(isempty(badMetrics), 'Metric not recognized: %s\n', badMetrics{:});
				else
					metrics = validMetrics(1:end-1);
				end
			end
		end
		
		
		function values = getValues(self, varargin) %% TODO ADD ABILITY TO SUB FIRST
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
				rawVals = reshape(unique(self.sampleMap.(param{:}), 'stable'), 1, []);
				if iscell(rawVals) && ischar(rawVals{1})
					% Get rid of empty strings
					values.(param{:}) = rawVals(~cellfun(@isempty, rawVals));
				else
					% Empty cells are not carried over for numeric values
					values.(param{:}) = rawVals;
				end
			end
		end
		
		
		function sampleIDs = getSampleIDs(self, treatments)
			% Returns an array of sample IDs corresponding with the given treatments 
			% in the order requested
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
			if (numel(numParams) > 1) && false
				sampleIDs = reshape(linearSampleIDs, numParams);
			else
				sampleIDs = linearSampleIDs;
			end
		end
		
		
		function [dataMatrix, gateLogicals] = slice(self, sampleIDs, sliceParams)
			% Slices the data, returning an N x M matrix of data from N cells in M channels. 
			%
			%	dataMatrix = self.slice(sampleID, sliceParams)
			%
			%	Inputs
			%		sampleIDs		<integer> The sample(s) to slice as given by
			%						the numerical sample ID(s).
			%
			%		sliceParams		<struct> Optional, struct with optional fields:
			%						'channels':  <cell, char>, defaults to self.channels
			%						'dataType':  <char>, defaults to 'raw'
			%						'gate':		 <char>, defaults to no gate
			%						'equalize':  <logical>, TRUE returns an equal
			%									 number of points from each sample
			%									 (default = FALSE), operation
			%									 performed before applying bins
			%						'numPoints': The number of cells to extract
			%									 per sample (if the minimum number 
			%									 of cells per sample is lower, the 
			%									 method will use that value). 
			%						'controls':  <logical> TRUE slices from
			%									 controlData rather than sampleData
			%									 (default = FALSE)
			%						'bins':		 <numeric> defaults to all cells
			%									 An Nx1 set of numerical bin IDs or an 
			%									 NxD set of bin coordinates where D = #
			%									 of bin channels. 
			%									 Automatically forces 'dataType' to be 
			%									 'self.binDataType' regardless of whether 
			%									 they are given or not
			%									 <Can input 'all' to select all bins>
			%
			%	Ouputs
			%		dataMatrix		<double> N x M matrix of data from the given
			%						sample where N is the number of cells in the
			%						returned data and M is the number of
			%						channels requested. 
			%
			%		gateLogicals	<cell> S x 1 cell array of N x 1 logical vector 
			%						indicating which points of sample ID S are in the 
			%						given gates/bins. 
			
			% Check and extract inputs
			[sliceData, sliceChannels, sliceDataType, sliceGates] = zCheckInputs_slice(self);
			
			% Slice out data
			dataMatrix = [];
			gateLogicals = cell(size(sliceGates));
			for s = 1:numel(sampleIDs)
				
				% Extract data
				sID = sampleIDs(s);
				dataS = zeros(numel(sliceGates{s}), numel(sliceChannels));
				for ch = 1:numel(sliceChannels)
					if ismember(sliceChannels, Gating.SCATTER_CHANNELS)
						sdt = 'raw';
					else
						sdt = sliceDataType;
					end
					dataS(:, ch) = sliceData(sID).(sliceChannels{ch}).(sdt)(sliceGates{s});
				end
				dataMatrix = [dataMatrix; dataS];
				
				gateLogicals{s} = false(sliceData(sID).nObs, 1);
				gateLogicals{s}(sliceGates{s}) = true;
			end
			
			
			 % --- Helper Functions --- %
			
			
			function [sliceData, sliceChannels, sliceDataType, sliceGates] = zCheckInputs_slice(self)
				% Validates slice properteis and that the sampleID is valid
				
				% Check and update slice parameters as needed
				if exist('sliceParams', 'var')
					validateattributes(sliceParams, {'struct'}, {}, mfilename, 'sliceParams', 2);
				else
					sliceParams = struct();
				end
				
				% Slices data from the given channels
				if isfield(sliceParams, 'channels')
					sliceChannels = sliceParams.channels;
					if ischar(sliceChannels), sliceChannels = {sliceChannels}; end % Force cell
					sliceChannels = reshape(sliceChannels, 1, []); % Force row vector
					
					badChannels = setdiff(sliceChannels, [self.channels, Gating.SCATTER_CHANNELS]);
					assert(isempty(badChannels), ...
							'Channel not allowed: %s\n', badChannels{:});
				else
					sliceChannels = self.channels; % Default is all channels
				end
				
				% Slices data from the given dataset
				if (isfield(sliceParams, 'controls') && sliceParams.controls)
					sliceData = catstruct(self.controlDataScatter, self.controlData);
					sliceControls = true;
				else
					sliceData = catstruct(self.sampleDataScatter, self.sampleData);
					sliceControls = false;
				end
				
				% Ensure sampleID is a valid integer
				validateattributes(sampleIDs, {'numeric'}, {'positive'}, ...
					mfilename, 'sampleIDs', 1);
				sampleIDs = unique(round(sampleIDs)); 
				assert(all(sampleIDs <= numel(sliceData)), ...
					'At least one sampleID is too large!')
				
				% Slices data of the given dataType
				if isfield(sliceParams, 'dataType')
					validatestring(sliceParams.dataType, self.dataTypes, mfilename, 'sliceParams.dataType');
					sliceDataType = sliceParams.dataType;
				else
					sliceDataType = 'raw'; % Default is raw data
				end
				
				% Applies the given gate to the sliced data
				if isfield(sliceParams, 'gate')
					validatestring(sliceParams.gate, self.gateNames, mfilename, 'sliceParams.gate');
					sliceGates = cell(1, numel(sampleIDs));
					for id = 1:numel(sampleIDs)
						sampleID = sampleIDs(id);
% 						sliceGates{id} = sliceData(sampleID).gates.(sliceParams.gate);
						sliceGates{id} = find(sliceData(sampleID).gates.(sliceParams.gate));
					end
				else
					for id = 1:numel(sampleIDs)
						sampleID = sampleIDs(id);
						numSliceCells = numel(sliceData(sampleID).(sliceChannels{1}).raw);
% 						sliceGates{id} = true(numSliceCells, 1);
						sliceGates{id} = (1:numSliceCells)'; % Default is all cells
					end
				end
				
				% Slice the same number of points from each sample
				numPoints = min(cellfun(@numel, sliceGates));
				subSample = false;
				if (isfield(sliceParams, 'numPoints') && ~isempty(sliceParams.numPoints))
					numPoints = min(numPoints, sliceParams.numPoints);
					subSample = true;
				end
				if ((isfield(sliceParams, 'equalize') && sliceParams.equalize) || subSample)
					% If sP.equalize TRUE and sP.numPoints not given, then use minPoints
					% If both or just sP.numPoints given, use min of numPoints sP.numPoints
					% If neither are given, then subsampling is completely skipped
					for id = 1:numel(sampleIDs)
						ss = FlowAnalysis.subSample(numel(sliceGates{id}), numPoints);
						sliceGates{id} = sliceGates{id}(ss);
					end
				end
				
				% Slices from the given bins
				if (isfield(sliceParams, 'bins') && ~isempty(self.bins))
					% Controls are not binned, so if they are the sliceData, we
					% just skip the rest of the inputs checking
					if sliceControls
						warning('Controls are not binned! Skipping')
						return
					end
					
					validateattributes(sliceParams.bins, {'numeric', 'char'}, {}, mfilename, 'sliceParams.bins');
					if ischar(sliceParams.bins)
						if strcmpi(sliceParams.bins, 'all')
							sliceParams.bins = (1:self.numBins)';
						else
							error('Bin char input not recognized: %s', sliceParams.bins)
						end
					end
					assert(size(sliceParams.bins, 2) == 1 || ...
						   size(sliceParams.bins, 2) == numel(fieldnames(self.binInputs)), ...
						   'Bin IDs formatted incorrectly!')
					assert(size(sliceParams.bins, 1) <= self.numBins, ...
						   'Too many bins requested!')
					
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
					if ~strcmpi(sliceDataType, self.binDataType)
						warning('Setting ''dataType'' to %s to match binning data', self.binDataType)
					end
					sliceDataType = self.binDataType;

					% Extract cells in the requested bins
					for id = 1:numel(sampleIDs)
						sampleID = sampleIDs(id);
						cellsInBins = [];
						for b = binIdxs % Sequentially looks at each bin ID
							cellsInBins = [cellsInBins, self.bins{sampleID}{b}];
						end

	% 						inBin = false(size(sliceGates{id}));  % Make 'gate' for cells in bins
	% 						inBin(cellsInBins) = true;		 % Fill out logical index array
	% 						sliceGates{id} = (sliceGates{id} & inBin); % Combine w/ sliceGate for simplicity
						sliceGates{id} = intersect(sliceGates{id}, cellsInBins);
					end
				end
			end
		end
		
		
		function self = operate(self, operation, values, initDataType, newDataType)
			% Applies the given operation to the given dataType, creating a new
			% dataType which is a mathematically transformed with the operation
			%
			%	FlowData.operate(operation, values, initDataType, newDataType)
			%
			%	Inputs
			%
			%		operation		<char> The mathematical operation to perform.
			%						  'add':	Adds the values to the data 
			%						  'scale':	Scales the data by the values
			%						  'exp':	Takes the data the value exponents
			%						  'power':	Takes the values to the data exponent
			%
			%		values			<numeric> A 1xC vector of values used for
			%						the operations, where C := # channels 
			%						* Can alternatively pass a single
			%						  value to apply to all channels
			%
			%		initDataType	<char> (Optional) The dataType to use for
			%						the operation (default = 'raw').
			%
			%		newDataType		<char> (Optional) The name of the new dataType 
			%						created by the operation (default = 'op#')
			%						where # indicates how many times operate has
			%						been called with no newDataType given. 
			
			zCheckInputs_operate(self)
			
			for iii = 1:max(self.numControls, self.numSamples)
				for ch = 1:numel(self.channels)
					
					chan = self.channels{ch};
					
					if (iii <= self.numControls && ~isempty(self.controlData(iii).(chan)))
						data = self.controlData(iii).(chan).(initDataType);
						self.controlData(iii).(chan).(newDataType) = applyOps(data, ch);
					end
					
					if iii <= self.numSamples
						data = self.sampleData(iii).(chan).(initDataType);
						self.sampleData(iii).(chan).(newDataType) = applyOps(data, ch);
					end
				end
			end
			
			self.addDataTypes(newDataType);
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_operate(self)
				
				validOps = {'add', 'scale', 'exp', 'power'};
				validatestring(operation, validOps, mfilename, 'operation', 1);
				validateattributes(values, {'numeric'}, {}, mfilename, 'scalar', 2);
				if (length(values) == 1)
					values = repmat(values, 1, size(self.channels));
				end
				assert(length(values) == numel(self.channels), ...
					'Incorrect scalar length passed!');
				
				if exist('initDataType', 'var')
					validatestring(initDataType, self.dataTypes, mfilename, 'initDataType', 3);
				else
					initDataType = 'raw';
				end
				
				if exist('newDataType', 'var')
					validateattributes(newDataType, {'char'}, {}, mfilename, 'initDataType', 3);
				else
					self.unnamedOps = self.unnamedOps + 1;
					newDataType = ['op', num2str(self.unnamedOps)];
				end
				
			end
			
			
			function opData = applyOps(data, ch)
				switch operation
					case 'add'
						opData = data + values(ch);
					case 'scale'
						opData = data .* values(ch);
					case 'exp'
						opData = data .^ values(ch);
					case 'power'
						opData = values(ch) .^ data;
				end
			end
			
		end
				
		
		function self = threshGate(self, channels, mode, thresh)
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
			
			zCheckInputs_threshGate(self);
			
			% Setup new directory for gates
			gateDir = [self.folder, 'Gating', filesep];
			if ~exist(gateDir, 'file')
				mkdir(gateDir)
			end
			gatesSaveName = [self.date, '_', self.name, '_Thresh_'];
						
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
					
					ss = FlowAnalysis.subSample(numel(combData), 1e4);
					
					figThresh = figure();
					ax = gca(); hold(ax, 'on')
					histogram(ax, Transforms.lin2logicle(combData(ss), false, self.logicleParams))
					title('Draw a line to set an x-axis threshold', 'fontsize', 16)
					ylabel('Count', 'fontsize', 14)
					xlabel(strrep(chan, '_', '-'), 'fontsize', 14)
					Plotting.biexpAxes(ax, true, false, false, false, self.logicleParams);
					
					h = imline();
					position = wait(h);
					
					saveas(figThresh, [gateDir, gatesSaveName, chan, '.fig']);
					
					% Fit the points to get a line, then use that to find the
					% x-intercept (threshold value)
					fit = polyfit(position(:, 1), position(:, 2), 1);
% 					plot(ax, (0:0.5:4.5), fit(1) * (0:0.5:4.5) + fit(2), 'r-');
					threshVal = Transforms.logicle2lin(-fit(2) / fit(1), false, self.logicleParams);
				end
				
				% Find cells passing threshold and record gate
				thrGateNames{ch} = ['TH', self.SHORT_COLORS.(chan)];
				for i = 1:self.numSamples
					passThresh = (self.sampleData(i).(chan).raw >= threshVal);
					self.sampleData(i).gates.(thrGateNames{ch}) = passThresh;
				end
				
				% Do the same for controls
				for i = 1:self.numControls
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
			
			
			function zCheckInputs_threshGate(self)
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
				
		
		function self = crossGates(self, gates, mode)
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
					for i = 1:self.numControls
						if ~isempty(self.controlData(i).(self.channels{1})) % Some tcData will be empty 
							self.controlData(i) = crossOR(self.controlData(i), newGateName);
						end
					end
				case 'and'
					for i = 1:self.numSamples
						self.sampleData(i) = crossAND(self.sampleData(i), newGateName);
					end
					for i = 1:self.numControls
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
			
			% Some gates may not be applied to the sample data, so only look at
			% names of gates literally assigned under sampleData
			sampleGateNames = fieldnames(self.sampleData(1).gates);
			numGates = numel(sampleGateNames);
			gatePcts = zeros(self.numSamples, numGates);
			
			for g = 1:numGates
				gate = sampleGateNames{g};
				
				for si = 1:self.numSamples
					gatePcts(si, g) = mean(self.sampleData(si).gates.(gate)) * 100;
				end
			end
			
			% Convert to readable table
			sampleGateNames = fieldnames(self.sampleData(1).gates);
			gatePcts = array2table(gatePcts);
			gatePcts.Properties.VariableNames = sampleGateNames;
			
			% Append to Sample Map
			gatePcts = [self.sampleMap, gatePcts];
			
			if ~ismember('noPrint', options)
				disp(gatePcts)
			end
		end
		
		
		function self = editSampleMap(self, newSampleMap)
			% Allows the user to change the sample map by passing a new one to
			% the object. This method ensures that the sample map is valid.
			%
			%	Inputs
			%		newSampleMap		<Table> A new sample map for the data
			
			validateattributes(newSampleMap, {'table', 'cell', 'char'}, {}, mfilename, 'newSampleMap', 1);
			
			% Handle non-table inputs
			if ischar(newSampleMap)
				if ~contains(newSampleMap, filesep)
					newSampleMap = [self.folder, newSampleMap];
				end
				assert(logical(exist(newSampleMap, 'file')), 'New Sample Map file not found!');
				newSampleMap = readtable(newSampleMap, 'delimiter', '\t');
			elseif iscell(newSampleMap)
				assert(self.sampleMap.width == size(newSampleMap, 2), ...
					'New Sample Map has incorrect number of variables (%d) compared to the data (%d)', ...
					size(newSampleMap, 2), self.sampleMap.width)
				varNames = self.sampleMap.Properties.VariableNames;
				newSampleMap = cell2table(newSampleMap);
				newSampleMap.Properties.VariableNames = varNames;
			end
			
			% Check table size
			assert(height(newSampleMap) == self.numSamples, ...
				'New Sample Map has incorrect number of samples (%d) compared to the data (%d)', ...
				height(newSampleMap), self.numSamples);
			
			% Re-assign variable
			self.sampleMap = newSampleMap;
			
			fprintf(1, 'Finished adding new sample map\n');
		end
		
		
		function self = editlogicleParams(self, newParams)
			% Allows the user to change the logicle conversion parameters used
			% by the object by passing a new set of parameters. This method
			% ensures that the new parameters are valid.
			
			validateattributs(newParams, {'struct'}, {}, mfilename, 'newParams', 1);
			
			% Look for any missing fields - keep existing parameters for those
			% that are not given
			if ~isfield(newParams, 'T'), newParams.T = self.logicleParams.T; end
			if ~isfield(newParams, 'M'), newParams.M = self.logicleParams.M; end
			if ~isfield(newParams, 'r'), newParams.r = self.logicleParams.r; end
			if ~isfield(newParams, 'MEF'), newParams.MEF = self.logicleParams.MEF; end
			
			% Check the values are valid
			validateattributes(newParam.T, {'numeric'}, {'scalar', 'positive'}, mfilename, 'newParams.T');
			validateattributes(newParam.M, {'numeric'}, {'scalar', 'positive'}, mfilename, 'newParams.M');
			validateattributes(newParam.r, {'numeric'}, {'scalar', 'negative'}, mfilename, 'newParams.r');
			validateattributes(newParam.MEF, {'numeric'}, {'scalar', 'positive'}, mfilename, 'newParams.MEF');
			
			self.logicleParams = newParams;
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
% 				case {'binInputs', 'binDataType'}
% 					self.bin(self.binInputs, self.binDataType)
% 				
% 				case {'test'}
% 					fprintf(1, 'Test set to: %s\n', event.AffectedObject.test);
% 			end
% 		end
		
		
		function self = addDataTypes(self, dataTypes)
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
		
		
		function self = addGates(self, gates)
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
		
		
		function filenames = convertToFullFile(~, filenames, folder)
			% Checks if a file is a full path and adds the path using the given
			% folder if not.
			%
			%	filenames: string or cell array of strings, folder: string
			
			% Convert to cell array for simplicity
			wasStr = false;
			if ~iscell(filenames), filenames = {filenames}; wasStr = true; end
			
			for fn = 1:numel(filenames)
				if ~strfind(filenames{fn}, filesep)
					filenames{fn} = [folder, filenames{fn}];
				end
			end
			
			% Convert back to string if that was provided
			if wasStr, filenames = filenames{1}; end
		end
	end
end