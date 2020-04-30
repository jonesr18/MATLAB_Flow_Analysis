classdef FlowData < matlab.mixin.Copyable
	% A data structure for managing flow cytometry data. Uses and references
	% several other classes and functions in the repository. 
	%
	%	Visible Properties
	%
	%		name			 <char>		Experiment name
	%		date			 <char>		Experiment start date
	%		folder			 <char>		The full path name for the experiment
	%		controlFolder	 <char>		The full path name for the controls
	%		cytometer		 <char>		The cytometer used for the experiment
	%		channels		 <cell>		Cytometer channel names
	%		colors			 <cell>		Fluorescent color names
	%
	%		sampleData		 <struct>	Sample fluorescence and gate data in standard struct
	%		controlData		 <struct>	Similar to sampleData but for controls.
	%									The order of controls should be the same
	%									as their corresponding channels. 
	%		numSamples		 <numeric>	The number of data samples
	%		numControls		 <numeric>	The number of control samples
	%		numCells		 <numeric>	The number of cells in each sample
	%		numCellsControls <numeric>	The number of cells in each control sample
	%		sampleMap		 <table>	Experimental information for samples
	%
	%		dataTypes		 <cell>		Cell array of data types 
	%									('raw', 'comp', 'mefl', etc)
	%		gateNames		 <cell>		Cell array of gate names (strings)
	%		gatePolygons	 <struct>	Mapping between gate names and polygons for sampleData
	%		threshGateVals	 <numeric>	Threshold gate thresh values
	%
	%		coefficients	 <numeric>	Compensation coefficients
	%		autfluor		 <numeric>  Autofluorescence per channel
	%		compDataType	 <char>		Data type used for compensation
	%
	%		mefFitsSamples	 <table>	Table containing MEF unit fits for samples
	%		mefFitsControls	 <table>	Table containing MEF unit fits for controls
	%		meflFits		 <struct>	Mapping between channel names to MEF-MEFL conversion factors 
	%		meflScalars		 <struct>	Mapping between channel names to raw-MEFL conversion
	%
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
	%		binNames		 <cell>		Named binning schemes for save/load
	%
	%	Public Methods
	%
	%		FlowData(dataFnames, expDetails)
	%		gate(self, options)
	%		customGate(self, gateName, sampleIDs, sliceParams, axProperties, gatePolygon)
	%		threshGate(self, channels, mode, thresh, options, sampleIDs)
	%		crossGates(self, gates, mode)
	%		printGatePcts(self, options)
	%		removeGates(self, gates)
	%		addControls(self, controlFolder, wildTypeFname, singleColorFnames, twoColorFnames)
	%		convertToMEF(self, beadsControls, beadsSamples, options)
	%		convertToMEFL(self, options)
	%		convertToMEFL(self, options)
	%		bin(self, binInputs, binDataType, sampleIDs, binName, binFuncs, doPar)
	%		loadBins(self, binName)
	%		computeStats(self, sampleIDs, sliceParams, metrics, options, minCells)
	%		getValues(self, varargin)
	%		getSampleIDs(self, treatments)
	%		slice(self, sampleIDs, sliceParams)
	%		operate(self, operation, values, initDataType, newDataType)
	%		editSampleMap(self, newSampleMap)
	%		editLogicleParams(self, newParams)
	%		rotate(self, channels, theta)
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
		controlFolder = '';			% Controls full path name
		cytometer = '';				% Experiment cytometer used
		channels = {};				% Cytometer channel names
		colors = {};				% Fluorescent color names
		
		sampleData = struct();		% Sample fluorescence and gate data in standard struct
		controlData = struct();		% Similar to sampleData but for controls
		numSamples = 0;				% The number of data samples
		numControls = 0;			% The number of controls
		numCells = [];				% The number of cells in each sample
		numCellsControls = [];		% The number of cells in each control
		sampleMap = table();		% Experimental information for samples
		
		dataTypes = {'raw'};		% Cell array of data types ('raw', 'mComp', 'mefl', etc)
		gateNames = {};				% Cell array of gate names (strings)
		gatePolygons = struct();	% Mapping between gate names and polygons for sampleData
		threshGateVals = [];		% Threshold gate thresh values
		
		coefficients = [];			% Compensation coefficients
		autofluor = [];				% Autofluorescence per channel
		compDataType = '';			% Data type used for compensation
		
		mefFitsSamples = table();	% Table containing MEF unit fits for samples
		mefFitsControls = table();	% Table containing MEF unit fits for controls
		meflFits = struct();		% Mapping between channel names to MEF-MEFL conversion factors 
		meflScalars = struct();		% Mapping between channel names to raw-MEFL conversion factors 
		
		bins = {};					% Cell array where each element corresponds with a data sample
		numBins = 0;				% Number of bins
		binSizes = [];				% Number of bins in each dimension
		binInputs = struct();		% Struct w/ bin channels as fields and edges as values
		binDataType = '';			% Binned dataType ('raw', 'mComp', 'mefl', etc)
		binNames = {};				% Named binning schemes for save/load
		
		logicleParams = struct( ... % Logicle transformation parameters to use
			'T', 2^18, ...
			'M', 4.5, ...
			'r', -150);
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
		
		function self = FlowData(dataFnames, expDetails)
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
			%		expDetails		<struct> A struct with the following fields
			%						recording experimental details:
			%							name		|	<char>
			%							date		|	<char>
			%							folder		|	<char> (full-path folder name)
			%						*	sampleMap	|	<char>
			%							cytometer	|	<char>
			%						**	channels	|	<char, cell>
			%						***	colors		|	<char, cell>
			%
			%						* 'sampleMap': name of the .txt file containing 
			%						treatments information for each sample. 
			%						 - See example file in source folder.
			%						** 'channels': The color channel(s) corresponding
			%						with the desired subset of data in dataStruct.
			%						 - Must be a field of dataStruct
			%						 - FSC/SSC are automatically taken
			%						*** 'colors': The name(s) of the fluorescent
			%						proteins or fluorophores that correspond
			%						with each measured channel
			%						 - Must match size and order of 'channels'
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
			self.channels = expDetails.channels;
			self.colors = expDetails.colors;
			
			% Extract sample map
			if ~isempty(sampleMapFname)
				try % Don't want to throw away loaded data
					sampleMap = readtable(sampleMapFname, 'delimiter', '\t');
					assert(height(sampleMap) == numel(dataStruct), ...
						'Sample map (%d) contains incorrect number of rows (%d)', ...
						height(sampleMap), numel(dataStruct))
					
					self.sampleMap = sampleMap;
				catch ME
					warning('Error occured while reading sample map, skipping...')
					disp(ME)
				end
			end
			
			% Extract data from the given channels
			self.numSamples = numel(dataStruct);
			self.numCells = zeros(self.numSamples, 1);
			for i = 1:self.numSamples
				
				nObs = dataStruct(i).nObs;
				self.numCells(i) = nObs;
				
				% Extract desired color channels
				for ch = [expDetails.channels, {'Time'}]
					self.sampleData(i).(ch{:}).raw = dataStruct(i).(ch{:}).raw;
				end
				self.sampleData(i).nObs = nObs;
				
				% Extract scatter channels
				for ch = Gating.SCATTER_CHANNELS
					self.sampleDataScatter(i).(ch{:}).raw = dataStruct(i).(ch{:}).raw;
				end
				self.sampleDataScatter(i).nObs = nObs;
			end
			
			% Add listeners for settable public properties
% 			addlistener(self, 'test', 'PostSet', @self.handlePropEvents);
% 			addlistener(self, {'binInputs', 'binDataType'}, ...
% 				'PostSet', @self.handlePropEvents);
			
			fprintf(1, 'Finished constructing FlowData object\n')
			
			
			% -- Helper functions -- %
			
			
			function [dataStruct, sampleMapFname] = zCheckInputs(self)
				validateattributes(dataFnames, {'cell', 'char'}, {}, mfilename, 'dataFilenames', 1);
				validateattributes(expDetails, {'struct'}, {}, mfilename, 'exptDetails', 2);
				
				% Check required experiment details are present
				requiredFields = {'date', 'name', 'folder', 'cytometer', 'channels', 'colors'};
				missingFields = setdiff(requiredFields, fieldnames(expDetails));
				assert(isempty(missingFields), 'Experiment details missing field: %s\n', missingFields{:});
				
				% Check channels is a valid input then force a cell array
				channels = expDetails.channels;
				validateattributes(channels, {'cell', 'char'}, {}, mfilename, 'expDetails.channels');
				if ischar(channels), channels = {channels}; end
				expDetails.channels = reshape(channels, 1, []); % Ensure row vector
				
				% Check colors is a valid input, force cell array, validate number
				colors = expDetails.colors;
				validateattributes(colors, {'cell', 'char'}, {}, mfilename, 'expDetails.channels');
				if ischar(colors), colors = {colors}; end
				expDetails.colors = reshape(colors, 1, []); % Ensure row vector
				
				% Verify all other expDetails fields are chars
				for fi = setdiff(fieldnames(expDetails), {'channels', 'colors'})'
					validateattributes(fi{:}, {'char'}, {}, mfilename, ['expDetails.', fi{:}])
				end
				
				% Check experiment folder exist and add filesep to end for later extensions
				assert(logical(exist(expDetails.folder, 'file')), 'Experiment folder does not exist');
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
				
				% Get data filenames and import Sample 1
				dataFnames = self.convertToFullFile(dataFnames, expDetails.folder);
				ds1 = FlowAnalysis.openFiles(dataFnames{1});
				
				% Check channels are present in dataStruct
				badChannels = setdiff(expDetails.channels, fieldnames(ds1));
				assert(isempty(badChannels), ...
					'Channel not in dataStruct: %s\n', badChannels{:});
				
				% Remove unwanted fields from structure
				goodFields = [expDetails.channels, Gating.SCATTER_CHANNELS, {'Time'}, {'nObs'}];
				ds1 = rmfield(ds1, setdiff(fieldnames(ds1), goodFields));
				
				% Import all data
				if (numel(dataFnames) > 1)
					ds = FlowAnalysis.openFiles(dataFnames{2:end});
					ds = rmfield(ds, setdiff(fieldnames(ds), goodFields));
					dataStruct = [ds1, ds];
				else
					dataStruct = ds1;
				end
			end
		end
		
		
		function self = gate(self, options)
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
			%
			%		options		<char, cell> (optional) Logical flags:
			%					'onlyP1': Only do P1 gating (FSC_A vs SSC_A)
			%					'recompute': Force a remake of the gates
						
			% Determine gating options
			onlyP1 = false;
			if exist('options', 'var')
				if islogical(options) % for reverse compatibility
					if options, options = {'onlyP1'}; else, options = {}; end
				elseif ischar(options)
					options = {options}; % Force to be a cell array
				end
				if ismember('onlyP1', options), onlyP1 = true; end
			else
				options = {};
			end
			self.onlyP1 = onlyP1; %#ok<*PROPLC>
			
			% Swap X/Y for LSRII-HTS2 since that one is always weird
			if strcmpi(self.cytometer, 'Koch-LSRII-HTS2')
				swap = true; 
			else
				swap = false; 
			end
			
			% Setup new directories for gates
			if (self.controlsAdded)
				gateDirControls = [self.controlFolder, 'Gating', filesep];
				if ~exist(gateDirControls, 'file')
					mkdir(gateDirControls)
				end
			end
			gateDirSamples = [self.folder, 'Gating', filesep];
			if ~exist(gateDirSamples, 'file')
				mkdir(gateDirSamples)
			end
			
			% Check if gates have already been made for this data first, then
			% process if necessary. 
			gatesSaveName = [self.date, '_', self.name];
			if (self.controlsAdded) 
				gatesFnameControls = [gateDirControls, 'Gate-Polygons_Controls.mat'];
				if (~ismember('recompute', options) && exist(gatesFnameControls, 'file'))
					% Load existing control gates
					load(gatesFnameControls, 'gateP1c', 'gateP2c', 'gateP3c');
				else
					% Do manual control gating
					[gateP1c, gateP2c, gateP3c, gateFigs] = Gating.standardGating(self.controlDataScatter, onlyP1, swap);
					save(gatesFnameControls, 'gateP1c', 'gateP2c', 'gateP3c');
					for f = fieldnames(gateFigs)'
						saveas(gateFigs.(f{:}), [gateDirControls, 'Gate-', f{:}, '_Controls']);
					end
				end
			end
			gatesFnameSamples = [gateDirSamples, 'Gate-Polygons_', gatesSaveName '.mat'];
			if (~ismember('recompute', options) && exist(gatesFnameSamples, 'file'))
				% Load existing sample gates
				load(gatesFnameSamples, 'gateP1s', 'gateP2s', 'gateP3s');
			else
				% Do manual sample gating
				[gateP1s, gateP2s, gateP3s, gateFigs] = Gating.standardGating(self.sampleDataScatter, onlyP1, swap);
				save(gatesFnameSamples, 'gateP1s', 'gateP2s', 'gateP3s');
				for f = fieldnames(gateFigs)'
					saveas(gateFigs.(f{:}), [gateDirSamples, 'Gate-', f{:}, '_', gatesSaveName]);
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
			if (self.controlsAdded)
				self.controlDataScatter = Gating.applyStandardGates(self.controlDataScatter, gateP1c, gateP2c, gateP3c, swap);
			end
			self.sampleDataScatter	= Gating.applyStandardGates(self.sampleDataScatter,  gateP1s, gateP2s, gateP3s, swap);
			
			% Transfer gate logicals to externally accesible data
			if (self.controlsAdded)
				for cd = 1:self.numControls
					self.controlData(cd).gates = self.controlDataScatter(cd).gates;
				end
			end
			for sd = 1:self.numSamples
				self.sampleData(sd).gates = self.sampleDataScatter(sd).gates;
			end
			
			fprintf(1, 'Finished standard gating\n')
		end
		
		
		function self = customGate(self, gateName, sampleIDs, sliceParams, axProperties, gatePolygon)
			% Creates a custom gate using data from the given samples 
			% and the given slice parameters.
			%
			%	self.customGate(gateName, sampleIDs, sliceParams)
			%
			%	Inputs
			%
			%		gateName		<char> The name of the new gate
			%							Avoid using 'x' as the first letter since
			%							we use 'xABC' to denote the opposite of gate 
			%							'ABC', though the code would still work. 
			%
			%		sampleIDs		<integer> (Optional) The sample(s) to slice as 
			%							given by the numerical sample ID(s).
			%								(Defaults to all cells)
			%
			%		sliceParams		<struct> (Optional) Struct with optional fields:
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
			%		axProperties	<struct> (Optional) Settable axes properties
			%
			%		gatePolygon		<numeric> (Optional) An array of gate polygon vertices 
			
			% Check inputs
			zCheckInputs_customGate(self)
			
			gateDirSamples = [self.folder, 'Gating', filesep];
			if ~exist(gateDirSamples, 'file')
				mkdir(gateDirSamples)
			end
			
			% Extract data and do gating
			if exist('gatePolygon', 'var')
				gateFig = [];
			else
				outData = self.slice(sampleIDs, sliceParams);
				[~, gatePolygon, gateFig] = Gating.gatePolygon( ...
						outData(:, 1), outData(:, 2), axProperties);
			end
			
			% Store gate info
			if (isfield(sliceParams, 'controls') && sliceParams.controls)
				for ci = sampleIDs
					if isempty(self.controlData(ci).(self.channels{1}))
						continue % Handle empty tcData
					end
					% Extract sample data and gate
					[outData, sliceGates] = self.slice(ci, sliceParams);
					inGate = Gating.gatePolygon(outData(:, 1), outData(:, 2), ...
							axProperties, gatePolygon);
					
					% Adjust index based on input gate
					inGateFixed = Gating.fixGateIdxs(sliceGates{1}, inGate);
					self.controlData(ci).gates.(gateName) = inGateFixed;
				end
			else
				for si = sampleIDs
					% Extract sample data and gate
					[outData, sliceGates] = self.slice(si, sliceParams);
					inGate = Gating.gatePolygon(outData(:, 1), outData(:, 2), ...
							axProperties, gatePolygon);
					
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
			
			
			function zCheckInputs_customGate(self)
				
				assert(ischar(gateName), 'Gate name must be a character vector');
				
				if exist('sampleIDs', 'var')
					validateattributes(sampleIDs, {'numeric'}, {'positive'}, mfilename, 'sampleIDs', 2);
					sampleIDs = reshape(sampleIDs, 1, []); % For-loop optimized
					if (isfield(sliceParams, 'controls') && sliceParams.controls)
						assert(max(sampleIDs) <= self.numControls, ...
								'Sample ID(s) are too large');
					else
						assert(max(sampleIDs) <= self.numSamples, ...
								'Sample ID(s) are too large');
					end
				elseif (isfield(sliceParams, 'controls') && sliceParams.controls)
					sampleIDs = 1:self.numControls;
				else
					sampleIDs = 1:self.numSamples;
				end
				
				if exist('sliceParams', 'var')
					validateattributes(sliceParams, {'struct'}, {}, mfilename, 'sliceParams', 3);
				else
					sliceParams = struct();
				end
				
				if exist('axProperties', 'var')
					validateattributes(axProperties, {'struct'}, {}, mfilename, 'axProperties', 4);
				else
					axProperties = struct();
				end
				
				if exist('gatePolygon', 'var')
					validateattributes(gatePolygon, {'numeric'}, {}, mfilename, 'gatePolygon', 5);
				end
			end
		end
		
		
		function self = removeGates(self, gates)
			% Removes the given gates if they exist in the data structure
			
			validateattributes(gates, {'cell', 'char'}, {}, mfilename, 'gates', 1);
			if ~iscell(gates), gates = {gates}; end % For simplicity
			
			% Add dataTypes
			removed = intersect(gates, self.gateNames);
			self.gateNames = setdiff(self.gateNames, removed, 'stable');
			
			% Remove gate data (works fine if 'removed' is empty)
			removePolygons = intersect(fieldnames(self.gatePolygons), removed);
			self.gatePolygons = rmfield(self.gatePolygons, removePolygons);
			for ci = 1:self.numControls
				if isempty(self.controlData(ci).(self.channels{1}))
					continue % Handle empty tcData
				end
				removeControls = intersect(fieldnames(self.controlData(ci).gates), removed);
				self.controlData(ci).gates = rmfield( ...
						self.controlData(ci).gates, removeControls);
			end
			for si = 1:numel(self.sampleData)
				removeSamples = intersect(fieldnames(self.sampleData(ci).gates), removed);
				self.sampleData(ci).gates = rmfield( ...
						self.sampleData(ci).gates, removeSamples);
			end
			
			if ~isempty(removed), fprintf(1, 'Removed gate: %s\n', removed{:}); end
		end
		
%{		
% 		function self = addControls2(self, controlFolder, sampleMap, colors)
% 			% TODO
% 			% Adds wild-type, single-color, and two-color (optional) data to the dataset
% 			% so that we can do compensation (single-colors) and MEFL conversion (two-colors).
% 			%
% 			% The method generates self.controlData, a struct array where single-color  
% 			% data from channel X is in position X, two-color data from channel X is
% 			% in position 2*X (if applicable), and wild-type data is in the last
% 			% position (regardless of the presence of two-color data). 
% 			% 
% 			%	self.addControls(controlFolder, wildTypeFname, singleColorFnames, twoColorFnames)
% 			%
% 			%	Inputs
% 			%		controlFolder	<char> The full-path folder name for controls
% 			%		
% 			%		sampleMap		<char> The filename which relates control
% 			%						file names to the fluorescent protein or
% 			%						fluorophore that they are a single- or two-
% 			%						color control for. (Including non-fluorescent 
% 			%						or wild-type cells)
% 			%		
% 			%		colors			<char, cell> The fluorescent protein or
% 			%						fluorophore names used in the experiment. 
% 			%						These must be represented in the sample map.
% 			
% 			[wildTypeData, singleColorData, twoColorData] = zCheckInputs_addControls(self);
% 			FITC_IDX = find(strcmpi('FITC_A', self.channels));
% 			self.controlFolder = controlFolder;
% 			
% 			% Extract data
% 			self.controlData = extractData([self.channels, {'nObs', 'Time'}], ...
% 						wildTypeData, singleColorData, twoColorData, FITC_IDX);
% 			
% 			% Extract scatter data
% 			self.controlDataScatter = extractData([Gating.SCATTER_CHANNELS, {'nObs'}], ...
% 						wildTypeData, singleColorData, twoColorData, FITC_IDX);
% 			
% 			% Find number of cells per control
% 			self.numControls = numel(self.controlData);
% 			self.numCellsControls = zeros(self.numControls, 1);
% 			for ci = 1:self.numControls
% 				nc = self.controlData.nObs;
% 				if isempty(nc), nc = 0; end
% 				self.numCellsControls = nc;
% 			end
% 			
% 			self.controlsAdded = true;
% 			fprintf(1, 'Finished adding controls\n');
% 			
% 			
% 			% --- Helper Functions --- %
% 			
% 			
% 			function [wildTypeData, singleColorData, twoColorData] = zCheckInputs_addControls(self)
% 				
% 				validateattributes(controlFolder, {'char'}, {}, mfilename, 'controlFolder', 1);
% 				assert(logical(exist(controlFolder, 'file')), 'Controls folder does not exist');
% 				if ~(controlFolder(end) == filesep), controlFolder = [controlFolder, filesep]; end
% 				
% 				validateattributes(wildTypeFname, {'cell', 'char'}, {}, mfilename, 'wildTypeFname', 2);
% 				validateattributes(singleColorFnames, {'cell', 'char'}, {}, mfilename, 'singleColorFnames', 3);
% 				
% 				% Convert to cell arrays if necessary for convenience
% 				if ischar(wildTypeFname), wildTypeFname = {wildTypeFname}; end
% 				if ischar(singleColorFnames), singleColorFnames = {singleColorFnames}; end
% 				
% 				% Check number of scFiles
% 				assert(numel(singleColorFnames) == numel(self.channels), ...
% 					'Incorrect number of single color controls');
% 				
% 				% Add full-path to filenames
% 				wildTypeFname = self.convertToFullFile(wildTypeFname, controlFolder);
% 				singleColorFnames = self.convertToFullFile(singleColorFnames, controlFolder);
% 				
% 				% Open files
% 				wildTypeData = FlowAnalysis.openFiles(wildTypeFname{:});
% 				singleColorData = FlowAnalysis.openFiles(singleColorFnames{:});
% 				
% 				% Add twoColorData if applicable
% 				if exist('twoColorFnames', 'var')
% 					validateattributes(twoColorFnames, {'cell', 'char'}, {}, mfilename, 'twoColorFnames', 4);
% 					if ischar(twoColorFnames), twoColorFnames = {twoColorFnames}; end
% 					assert(numel(twoColorFnames) == sum(~strcmpi('FITC_A', self.channels)), ...
% 						'Incorrect number of two color controls');
% 					twoColorFnames = self.convertToFullFile(twoColorFnames, controlFolder);
% 					twoColorData = FlowAnalysis.openFiles(twoColorFnames{:});
% 				else
% 					twoColorData = [];
% 				end
% 			end
% 			
% 			
% 			function outData = extractData(channels, wtData, scData, tcData, FITC_IDX)
% 				% Generates the new controlData struct based on the given
% 				% individual structs and the desired channels
% 				%
% 				%	The FITC_IDX input is needed to tell the function which
% 				%	channel ID to skip when extracting two-color data. There is
% 				%	no two-color data for FITC since the FITC channel itself is
% 				%	the reference for the other fluorescent proteins. 
% 				
% 				for ch = channels
% 					for sc = 1:numel(scData)
% 						outData(sc).(ch{:}) = scData(sc).(ch{:});
% 					end
% 					
% 					tc = 0;
% 					if ~isempty(tcData)
% 						tcIdx = 0;
% 						for tc = 1:numel(scData) % Go over length of scData since they should match
% 							if tc == FITC_IDX
% 								% No two color controls for FITC channel, since it
% 								% is the reference color for MEFL conversion
% 								outData(sc + tc).(ch{:}) = [];
% 							else
% 								tcIdx = tcIdx + 1; % For indexing tcData separate of tc iterator
% 								outData(sc + tc).(ch{:}) = tcData(tcIdx).(ch{:});
% 							end
% 						end
% 					end
% 					outData(sc + tc + 1).(ch{:}) = wtData.(ch{:});
% 				end
% 			end
% 		end
%}
		
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
			%			%
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
			self.controlData = extractData([self.channels, {'Time'}], ...
						wildTypeData, singleColorData, twoColorData, FITC_IDX);
			
			% Extract scatter data
			self.controlDataScatter = extractData(Gating.SCATTER_CHANNELS, ...
						wildTypeData, singleColorData, twoColorData, FITC_IDX);
			
			% Find number of cells per control
			self.numControls = numel(self.controlData);
			self.numCellsControls = zeros(self.numControls, 1);
			for ci = 1:self.numControls
				nc = self.controlData(ci).nObs;
				if isempty(nc), nc = 0; end
				self.numCellsControls(ci) = nc;
			end
			
			self.controlsAdded = true;
			fprintf(1, 'Finished adding controls\n');
			
			
			% --- Helper Functions --- %
			
			
			function [wildTypeData, singleColorData, twoColorData] = zCheckInputs_addControls(self)
				
				validateattributes(controlFolder, {'char'}, {}, mfilename, 'controlFolder', 1);
				assert(logical(exist(controlFolder, 'file')), 'Controls folder does not exist');
				if ~(controlFolder(end) == filesep), controlFolder = [controlFolder, filesep]; end
								
				validateattributes(wildTypeFname, {'cell', 'char'}, {}, mfilename, 'wildTypeFname', 2);
				validateattributes(singleColorFnames, {'cell', 'char'}, {}, mfilename, 'singleColorFnames', 3);
				
				% Convert to cell arrays if necessary for convenience
				if ischar(wildTypeFname), wildTypeFname = {wildTypeFname}; end
				if ischar(singleColorFnames), singleColorFnames = {singleColorFnames}; end
				
				% Check number of scFiles
				assert(numel(singleColorFnames) == numel(self.channels), ...
						'Incorrect number of single color controls');
				
				% Add full-path to filenames if not included
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
						if ~isfield(scData(sc), ch{:})
							warning('Fieldname %s not found in sc control file %d', ch{:}, sc)
						else
							outData(sc).(ch{:}).raw = scData(sc).(ch{:}).raw;
						end
						
						if ~isfield(outData(sc), 'nObs')
							outData(sc).nObs = scData(sc).nObs;
						end
					end
					
					tc = 0;
					if ~isempty(tcData)
						tcIdx = 0;
						for tc = 1:numel(scData) % Go over length of scData since they should match
							if tc == FITC_IDX
								% No two color controls for FITC channel, since it
								% is the reference color for MEFL conversion
								outData(sc + tc).(ch{:}) = [];
								outData(sc + tc).nObs = [];
							else
								tcIdx = tcIdx + 1; % For indexing tcData separate of tc iterator
								outData(sc + tc).(ch{:}).raw = tcData(tcIdx).(ch{:}).raw;
								outData(sc + tc).nObs = tcData(tcIdx).nObs;
							end
						end
					end
					outData(sc + tc + 1).(ch{:}).raw = wtData.(ch{:}).raw;
					outData(sc + tc + 1).nObs = wtData.nObs;
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
			%	Note that all channels initially use the MEFL bead counts for
			%	fitting (since not all channels have bead units and the
			%	relationship between equivalent fluorophore counts are linear).
			%	These MEFL units are _unadjusted_ and need to be MEFL converted
			%	after running this method to 
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
			%			recompute		Flag to force recalculation
			
			assert(self.controlsAdded, 'Controls must be added before converting to MEF units\n');
			
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
			colorChans = reshape([self.colors; self.channels], 1, []);
			colorChans = cellfun(@(x) strrep(x, '_', '-'), colorChans, 'uniformoutput', false);
			ccString = sprintf('_%s', colorChans{:});
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
				mefFname = [beadDirs{b}, 'MEF-Fits_', beadSaveName, ccString, '.mat'];
				
				if ~ismember('recompute', options) && logical(exist(mefFname, 'file'))
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
									[beadDirs{b}, 'MEF-Fit_', self.colors{chID}, '_', ...
									strrep(self.channels{chID}, '_', '-'), '_', ...
									MEF_units{chID}, '_', beadSaveName, '.fig']); 
						end
						if isfield(figFits, 'manualPeaks')
							saveas(figFits.manualPeaks, [beadDirs{b}, 'Manual-Peaks_', beadSaveName, '.fig']);
						end
						if isfield(figFits, 'gate')
							if (b == 1) % Control beads
								gateDir = [self.controlFolder, 'Gating', filesep];
							else % Sample beads
								gateDir = [self.folder, 'Gating', filesep];
							end
							saveas(figFits.gate, [gateDir, 'Gate-Beads_', beadSaveName, '.fig']);
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
			self.mefFitsControls = fitsControls;
			self.mefFitsSamples = fitsSamples;
			
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
					'Control beads cytometer does not match experiment cytometer')
				assert(strcmpi(beadsSamples.cytometer, self.cytometer), ...
					'Sample beads cytometer does not match experiment cytometer')
				
				if ~exist('options', 'var')
					options = {};
				end
			end
		end
		
		
		function self = convertToMEFL(self, options)
			% Converts each channel to MEFL units using the compensated MEF units. 
			% The two-color controls are utilized to get ratios between each MEF
			% unit and MEFLs. Conversion factors are stored in self.meflConversions.
			%
			%	self.convertToMEFL(showPlots);
			%
			%	Adds new dataTypes: {'mefl'}
			%
			%	Inputs
			%		options			<cell> Contains optional string flags:
			%			showPlots		Flag to show fitted plots
			%			recompute		Flag to force recalculation

			% Check pre-requisites for running
			assert(self.controlsAdded, 'Controls must be added before converting to MEFL units\n');
			assert(self.mefConverted, 'MEF conversion must be run before converting to MEFL units\n');
			
			zCheckInputs_convertToMEFL();
			
			% Check if conversions already exist
			beadDir = [self.controlFolder, 'Calibration', filesep];
			if ~exist(beadDir, 'file')
				error('Controls bead directory not found. It should be set up during MEF calibration')
			end
			colors_channels = reshape([self.colors; strrep(self.channels, '_', '-')], 1, []);
			colorString = sprintf('_%s', colors_channels{:});
			meflFname = [beadDir, 'MEFL-Fits', colorString, '.mat'];
			
			if ~ismember('recompute', options) && logical(exist(meflFname, 'file'))
				fprintf(1, 'Loading pre-computed MEFL conversions\n');
				load(meflFname, 'meflFitsCalc')
			else
				% Compute mefl conversions
				tcData = self.controlData(numel(self.channels) + 1 : 2 * numel(self.channels));
				[meflFitsCalc, figFits] = Transforms.calibrateMEFL(tcData, self.channels, ...
						self.colors, 'mef', ismember('showPlots', options));
			
				save(meflFname, 'meflFitsCalc');
				
				FITC_IDX = find(strcmpi('FITC_A', self.channels));
				if ~isempty(fieldnames(figFits))
					for chID = setdiff(1:numel(self.channels), FITC_IDX)
						saveas(figFits.(self.colors{chID}), ...
								[beadDir, 'MEFL-Fit_', self.colors{chID}, '_', ...
								strrep(self.channels{chID}, '_', '-'), '_', ...
								self.colors{FITC_IDX}, '_', ...
								strrep(self.channels{FITC_IDX}, '_', '-'), '.fig']); 
					end
				end
			end
			
			% Add converted data as 'mefl' data type
			scalars = struct;
			for ci = 1:numel(self.channels)
				
				chan = self.channels{ci};
				color = self.colors{ci};
				
				% Add MEFLs for controls
				for si = 1:self.numControls
					if isempty(self.controlData(si).(chan)), continue, end % Some tcData will be empty 
					self.controlData(si).(chan).mefl = self.controlData(si).(chan).mef * meflFitsCalc.(color);
				end
				
				% Add MEFLs for sample
				for si = 1:self.numSamples
					self.sampleData(si).(chan).mefl = self.sampleData(si).(chan).mef * meflFitsCalc.(color);
				end
				
				scalars.(chan) = meflFitsCalc.(color) * 10^self.mefFitsSamples.(chan)(2);
			end
			
			self.meflFits = meflFitsCalc; 
			self.meflScalars = scalars;
			self.addDataTypes('mefl');
			self.meflConverted = true;
			fprintf(1, 'Finished converting to MEFL\n');
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_convertToMEFL()
				if ~exist('options', 'var')
					options = {};
				else
					if ~iscell(options)
						options = {options};
					end
				end
			end
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
			colors_channels = reshape([self.colors; strrep(self.channels, '_', '-')], 1, []);
			colorString = sprintf('_%s', colors_channels{:});
			compFname = [compDir, 'Comp-Fits_', dataType, colorString, '.mat'];
			if (exist(compFname, 'file') && ~all(logical(options.recompute)))
				% Load existing sample gates
				fprintf(1, 'Loading pre-computed coefficients and autofluorescence\n');
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
						scData, self.channels, self.colors, options); 
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
						scData, self.channels, self.colors, options);
			end
			
			% Save data/figures (default = do saving)
			if (~isfield(options, 'save') || options.save)
				save(compFname, 'coeffs', 'ints')
				if (isfield(fitFigs, 'pre') && ~isempty(fitFigs.pre)) 
					% Doesn't run if data was re-loaded or plots were not requested
					saveas(fitFigs.pre, [compDir, 'Pre-Comp_', dataType, colorString, '.fig'])
					saveas(fitFigs.post, [compDir, 'Post-Comp_', dataType, colorString, '.fig'])
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
				if ismember(dataType, {'mef', 'mefl'})
					options.scale = 1e3;
				else
					options.scale = 1;
				end
			end
			
		end
		
		
		function self = bin(self, binInputs, binDataType, sampleIDs, binName, binFuncs, doPar)
			% Sorts the sample data using the given set of channels/edges into a
			% number of bins using the given dataType for assignments. 
			%
			%	self.bin(binInputs, binDataType, sampleIDs, binName, binFuncs, doPar)
			%
			%	Inputs
			%		binInputs		<struct> A struct with channel names as keys and 
			%						bin edges as values. The struct tells the function 
			%						which channels to bin on and where to draw the bins 
			%						in each dimension. 
			%						To create bins with edges defined by the ratio 
			%						between two channels, create a substruct within 
			%						binInputs with 'ratio' in the key name and the
			%						following interface:
			%							'Channel1': channel_1_name (char)
			%							'Channel2': channel_2_name (char)
			%							'edges':	edges of ch1:ch2 ratios
			%						To create bins with edges defined by the product 
			%						between two channels, create a substruct within 
			%						binInputs with 'product' in the key name and the
			%						following interface:
			%							'Channel1': channel_1_name (char)
			%							'Channel2': channel_2_name (char)
			%							'edges':	edges of ch1:ch2 products
			%						NOTE: Cells with negative values in any channel
			%						are auto-discluded from ratio/product bins.
			%		
			%		binDataType		<char> The cell dataType to use (eg 'mefl', 'comp')
			%
			%		sampleIDs		<numeric> (Optional) This allows for constraining 
			%						which samples to bin (possibly saving time). 
			%						 - The default behavior is to bin all samples. 
			%						 - Pass an empty array to skip. 
			%
			%		binName			<char> (Optional) A name with which to identify
			%						the given binning scheme via saved files. If
			%						no name is given, defaults to a serial number 
			%						representing the Nth binning operation run.
			%
			%		binFuncs		<struct> (Optional) A struct with channel names as
			%						keys and edge functions as values. The channel names
			%						must match those in 'binInputs'. The functions are
			%						anonymous functions that allow altering the bin edges
			%						for a channel as a function of the other channels.
			%							Example: 
			%								struct('PE_A', @(x) x / (x(1) + 1e5))
			%							Here, x is assumed to be an array representing the
			%							channels passed into binInputs, in the order that
			%							they were entered into the struct (ie 
			%							fieldnames(binInputs)). So for two channels,
			%							'Pacific_Blue_A', and 'PE_A', the above code
			%							would adjust PE bins based on Pac Blue values. 
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
			
			[binFieldnames, binChannels, binEdges] = zCheckInputs_bin(self);
			
			slicedData = cell(1, self.numSamples);
			for sid = sampleIDs
				
				slicedData{sid} = self.slice(sid, struct( ...
						'channels', {binChannels}, ...
						'dataType', binDataType, ...
						'equalize', false));
				
			end
			
			binData = cell(size(slicedData));
			for sid = sampleIDs
				
				binData{sid} = zeros(size(slicedData{sid}, 1), numel(binFieldnames));
				
				for fbi = 1:numel(binFieldnames)
					
					fieldname = binFieldnames{fbi};
					
					if contains(fieldname, 'ratio', 'ignorecase', true)
						ratioFields = fieldnames(binInputs.(fieldname));
						ratioFieldsCh = ratioFields(contains(ratioFields, 'channel', 'ignorecase', true));
						
						% Initiate data as a ones vector to then multiply/divide
						% the ensuing input channels
						binData{sid}(:, fbi) = ones(size(slicedData{sid}, 1), 1);
						
						for rfci = 1:numel(ratioFieldsCh)
							chanName = ratioFieldsCh{rfci};
							chi = strcmpi(binChannels, binInputs.(fieldname).(chanName));
							if all(lower(chanName(1:4)) == 'inv_') % Inverted channel
								chanData = 1 ./ slicedData{sid}(:, chi);
							else
								chanData = 1 .* slicedData{sid}(:, chi);
							end
							binData{sid}(:, fbi) = binData{sid}(:, fbi) .* chanData;
						end
					else
						chi = strcmpi(binChannels, fieldname);
						binData{sid}(:, fbi) = slicedData{sid}(:, chi);
					end
				end
			end
			
			bins = cell(size(slicedData));
			if doPar
				parfor sid = sampleIDs
					bins{sid} = FlowAnalysis.simpleBin(binData{sid}, binEdges, binFuncs);
				end
			else
				for sid = sampleIDs
					bins{sid} = FlowAnalysis.simpleBin(binData{sid}, binEdges, binFuncs);
				end
			end
			
			% Extract size of each dimension independently, otherwise all
			% singular dimensions over 2D are lost! Also handles single-channel 
			% binning so only one size is given (rather than two).
			binSizes = zeros(1, numel(binChannels));
			for ci = 1:numel(binChannels)
				binSizes(ci) = size(bins{sampleIDs(1)}, ci);
			end
			
			% Only replace data for samples binned
			for sid = sampleIDs
				self.bins{sid} = bins{sid};
			end
			numBins = numel(bins{sampleIDs(1)});
			self.numBins = numBins;
			self.binSizes = binSizes;
			self.binInputs = binInputs;
			self.binDataType = binDataType;
			self.binned = true;
			self.binNames = [self.binNames, {binName}];
			
			% Setup new directory for bins
			binDir = [self.folder, 'Binning', filesep];
			if ~exist(binDir, 'file')
				mkdir(binDir)
			end
			binSaveName = [self.date, '_', self.name];
			binFname = [binDir, 'Binning_', binName, '_', binSaveName, '.mat'];
			
			save(binFname, 'bins', 'numBins', 'binSizes', 'binInputs', 'binDataType');
			
			
			fprintf(1, 'Finished binning\n')
			
			
			% --- Helper Functions --- %
			
			
			function [binFieldnames, binChannels, binEdges] = zCheckInputs_bin(self)
				% Validates that the given bin properties are ok, then sets the
				% object's properties themselves if all are ok.

				% Check properties
				validateattributes(binInputs, {'struct'}, {}, mfilename, 'inputs', 1);
				binFieldnames = fieldnames(binInputs)';
				binChannels = {};
				binEdges = cell(1, numel(binFieldnames));
				for fi = 1:numel(binFieldnames)
					fn = binFieldnames{fi};
					if contains(fn, 'ratio', 'ignorecase', true)
						rfn = fieldnames(binInputs.(fn));
						rfCh = rfn(contains(rfn, 'channel', 'ignorecase', true));
						for rci = 1:numel(rfCh)
							binChannels = [binChannels, {binInputs.(fn).(rfCh{rci})}];
						end
						binEdges{fi} = binInputs.(fn).edges;
					else % Normal channel input
						binChannels = [binChannels, {fn}];
						binEdges{fi} = binInputs.(fn);
					end
				end
				binChannels = unique(binChannels, 'stable');
				badChannels = setdiff(binChannels, self.channels);
				assert(isempty(badChannels), ...
						'Channel not allowed: %s\n', badChannels{:});
				
				validateattributes(binDataType, {'char'}, {}, mfilename, 'binDataType', 2);
				assert(any(strcmp(binDataType, self.dataTypes)), ...
						'Bin data type does not match any existing data types: %s\n', binDataType);
				
				if exist('sampleIDs', 'var') && ~isempty(sampleIDs)
					validateattributes(sampleIDs, {'numeric'}, {'positive'}, ...
							mfilename, 'sampleIDs', 3);
					sampleIDs = reshape(unique(round(sampleIDs)), 1, []); 
					assert(all(sampleIDs <= self.numSamples), ...
							'At least one sampleID is too large')
				else
					sampleIDs = 1:self.numSamples;
				end
				
				if exist('binName', 'var')
					validateattributes(binName, {'char'}, {}, mfilename, 'binName', 4);
				else
					binName = num2str(numel(self.binNames) + 1);
				end
				
				% Check binFuncs if applicable
				if exist('binFuncs', 'var')
					validateattributes(binFuncs, {'function_handle', 'cell'}, {}, mfilename, 'binFuncs', 5);
					% The rest of the validation occurs within FlowAnalysis.simpleBin()
				else
					% In the case of no binFuncs, all are independent
					binFuncs = cell(size(binEdges));
					binFuncs(:) = {@(x) 1};
				end
				
				doPar = exist('doPar', 'var') && all(logical(doPar));
			end
		end
		
		
		function self = loadBins(self, binName)
			% Loads bin IDs and other binning information. Defaults by drawing
			% from the auto-generated binning file, but can be given a specific
			% file to draw from as desired.
			%
			%	Inputs
			%
			%		binFname		<char> (Optional) A binName corresponding to
			%						a binning scheme previously run (see bin).
			%						If no input is given, the most recent binning
			%						scheme is loaded.
			
			zCheckInputs_loadBins(self);
			
			binDir = [self.folder, 'Binning', filesep];
			if ~exist(binDir, 'file')
				error('Bin dir not found')
			end
			binSaveName = [self.date, '_', self.name];
			binFname = [binDir, 'Binning_', binName, '_', binSaveName, '.mat'];
			
			if exist(binFname, 'file')
				ld = load(binFname);
			else
				error('Bin file not found')
			end
			
			if ismember('bins', fieldnames(ld))
				self.bins = ld.bins;
			else
				error('Missing variable: ''bins''.')
			end
			
			if ismember('binInputs', fieldnames(ld))
				self.binInputs = ld.binInputs;
			else
				error('Missing variable: ''binInputs''.')
			end
			
			if ismember('binDataType', fieldnames(ld))
				self.binDataType = ld.binDataType; 
			else
				error('Missing variable: ''binDataType''.')
			end
			
			sampleID = find(~isempty(self.bins), 1); % In case not all samples were binned
			self.binSizes = size(self.bins{sampleID});
			self.numBins = numel(self.bins{sampleID});
			
			fprintf(1, 'Finished loading binning data\n');
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_loadBins(self)			
				if ~exist('binName', 'var')
					if isempty(self.binNames)
						error('No binning schemes found - please provide a name')
					else
						binName = self.binNames{end};
					end
				end
			end
		end
		
		
		function stats = computeStats(self, sampleIDs, sliceParams, metrics, options, minCells)
			% Computes the given statistic metrics on the combined data 
			% from the given sampleIDs. Can get stats for independent bins. 
			%
			%	stats = self.computeStats(sampleIDs, sliceParams, metrics)
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
			%						'bins':		<numeric, char> One of several inputs:
			%									'none': Samples not split into bins
			%									'all': Samples divided into all bins
			%									<numeric>: An Nx1 set of numerical 
			%										bin IDs or an NxD set of bin 
			%										coordinates (D := bin channels) 
			%									If binning, automatically forces 
			%									'dataType' to be 'self.binDataType' 
			%									** Requesting anything but all bins
			%									will cause the output to be 2D rather 
			%									than N+1 D as the data will not be 
			%									reshaped to match the bin sizes.  
			%
			%		metrics			<char, cell> (Optional) A list of metrics to compute
			%							'numcells', 'pctpos', %iles: 'pX.Y' --> X.Y%, 
			%							'median', 'mean', 'geomean', 'stdev', 
			%							'geostdev', 'cv', 'sem', 'semb'*
			%							 * semb = bootstrapped SEM - slow!
			%							 -> If 'all' or no metrics are given, then
			%								all metrics except 'semb' are computed,
			%								along with %iles p10, p50, and p90. 
			%							 -> For %iles, the output names are 'pX_Y' 
			%								since '.' cannot be used in field names
			%
			%		options			<char, cell> (Optional) Optional inputs:
			%							'pos': Calculate metrics for only the 
			%							positive data - eg, useful for looking at
			%							smaller percentiles to plot on log-log.
			%							
			%		minCells		<numeric> (Optional) The minimum number of cells
			%							to calculate statistics (useful for binning).
			%
			%	Outputs
			%		stats		<struct> A struct where each field is the name of 
			%					a statistical metric and the value is a (N+1)D
			%					matrix of the metric values in N bin dimensions 
			%					across C channels (given in sliceParams) in the 
			%					final dimension. If no binning has been done,
			%					then the first dimension is singular. 
			%					--> In the case of requesting stats from specific bins
			%						rather than all of them, we do not structure the 
			%						bin stats into their "true" shape and instead return 
			%						a BxC matrix of metric values in each of the 
			%						requested B bins in all the C channels. 
			%
			%	TODO: Change output to matrix w/ last dimension metrics
			
			zCheckInputs_computeStats(self);
			
			% Set up binStats struct | Check size w/ 1st dim because that = # bins
			doBinning = ~isempty(sliceParams.bins); % Need since we overwrite sliceParams.bins
			if doBinning % Separate samples into bins
				sliceBins = sliceParams.bins;
				statSize = size(sliceBins, 1); 
				iterMax = size(sliceBins, 1);
			else
				statSize = 1;
				iterMax = 1;
			end
			
			% Set up output matrix
			for m = metrics
				stats.(m{:}) = nan([statSize, numel(sliceParams.channels)]);
			end
			stats.numcells = zeros([statSize, 1]);
			
			% Iterate over bins, compute stats for each one independently
			for i = 1:iterMax
				
				% Set bin to extract from
				if doBinning
					sliceParams.bins = sliceBins(i, :);
				end
				
				% Slice out data to compute stats with
				slicedData = self.slice(sampleIDs, sliceParams);
				stats.numcells(i) = size(slicedData, 1);
% 				size(slicedData, 1)
				
				if isempty(slicedData), slicedData = nan; end
				
				for ch = 1:numel(sliceParams.channels)
					
					% Compute linear index
					statIdx = i + iterMax * (ch - 1);
					
					% Extract channel data
					dataOut = slicedData(:, ch);
					dataOut = dataOut(~isnan(dataOut));
					posData = dataOut(dataOut > 0);
					if ismember('pctpos', metrics)
						chi = strcmpi(self.channels, sliceParams.channels{ch});
						thresh = self.threshGateVals(chi);
						stats.pctpos(statIdx) = 100 * sum(dataOut > thresh) ./ numel(dataOut);
					end
					if ismember('pos', options), dataOut = posData;	end
					if (numel(dataOut) < minCells), continue; end
					
					% Calculate percentiles
					pcts = setdiff(metrics(contains(metrics, 'p')), 'pctpos');
					if ~isempty(pcts)
						pctsNum = cellfun(@(x) str2double(x(2:end)), pcts);
						prctiles = prctile(dataOut, pctsNum);
						for pi = 1:numel(pcts)
							stats.(strrep(pcts{pi}, '.', '_'))(statIdx) = prctiles(pi);
						end
					end
					
					% Calculate other metrics
					if ismember('median', metrics), stats.median(statIdx) = median(dataOut); end % Same as p50, but preserving input name
					if ismember('mean', metrics), stats.mean(statIdx) = mean(dataOut); end
					if ismember('geomean', metrics), stats.geomean(statIdx) = geomean(posData); end
					if ismember('stdev', metrics), stats.stdev(statIdx) = std(dataOut, 0); end
					if ismember('geostdev', metrics), stats.geostdev(statIdx) = geostd(posData, 0); end
					if ismember('cv', metrics), stats.cv(statIdx) = mean(dataOut) / std(dataOut, 0); end
					if ismember('sem', metrics), stats.sem(statIdx) = std(dataOut, 0) / sqrt(numel(dataOut)); end
					if ismember('semb', metrics), stats.semb(statIdx) = semBootstrap(dataOut); end % NOTE: SUPER SLOW
% 					if ismember('ci95', metrics), binStats.CI95(binStatIdx) = ci95(dataOut); end
				end
			end
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_computeStats(self)
				
				validateattributes(sampleIDs, {'numeric'}, {'positive'}, mfilename, 'sampleIDs', 1);
				sampleIDs = reshape(unique(ceil(sampleIDs)), 1, []);
				assert(max(sampleIDs) <= self.numSamples, 'Sample ID(s) too large')
				
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
				
				% If no information about bins are given, do not separate
				% samples into bins. Otherwise, check if binning has been done
				% then figure out which bins to use. 
				if (~isfield(sliceParams, 'bins') || isempty(self.bins))
					sliceParams.bins = [];
				end
				if ischar(sliceParams.bins)
					if strcmpi(sliceParams.bins, 'all')
						sliceParams.bins = (1:self.numBins)';
					else
						sliceParams.bins = [];
					end
				end
				if ~isempty(sliceParams.bins)
					assert(size(sliceParams.bins, 2) == 1 || ...
						   size(sliceParams.bins, 2) == numel(fieldnames(self.binInputs)), ...
						   'Bin IDs formatted incorrectly')
					assert(size(sliceParams.bins, 1) <= self.numBins, ...
						   'Too many bins requested')
				end
				
				defaultMetrics = {'numcells', 'p10', 'p50', 'p90', 'mean', ...
						'geomean', 'stdev', 'geostdev', 'cv', 'sem'};
				if exist('metrics', 'var')
					validateattributes(metrics, {'cell', 'char'}, {}, mfilename, 'metrics', 3);
					if ischar(metrics), metrics = {metrics}; end % For simplicity
					metrics = reshape(metrics, 1, []); % Force row vector
					metrics = cellfun(@lower, metrics, 'uniformoutput', false);
				else
					metrics = defaultMetrics;
				end
				
				if exist('options', 'var')
					validateattributes(options, {'cell', 'char'}, {}, mfilename, 'options', 4);
					if ischar(options), options = {options}; end % For simplicity
					options = reshape(options, 1, []); % Force row vector
					options = cellfun(@lower, options, 'uniformoutput', false);
				else
					options = {};
				end
				
				if exist('minCells', 'var')
					validateattributes(minCells, {'numeric'}, {'scalar'}, mfilename, 'numCells', 5);
				else
					minCells = 1;
				end
			end
		end
		
		
		function values = getValues(self, varargin) %% TODO ADD ABILITY TO SUB FIRST
			% Returns unique values for each given experimental parameter.
			%
			%	values = self.getValues(parameters)
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
			if isempty(varargin)
				parameters = self.sampleMap.Properties.VariableNames;
			else
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
			end
			
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
			if (numel(numParams) > 1) && false % Currently not done
				sampleIDs = reshape(linearSampleIDs, numParams);
			else
				sampleIDs = reshape(linearSampleIDs, 1, []); % Force row vector for easier for-loop access
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
			%						'bins':		 <numeric> Defaults to all cells
			%									 An Nx1 set of numerical bin IDs or an 
			%									 NxD set of bin coordinates where D = #
			%									 of bin channels. 
			%									 Automatically forces 'dataType' to be 
			%									 'self.binDataType'
			%									 <Can input 'all' to select all bins>
			%
			%	Ouputs
			%		dataMatrix		<double> N x M matrix of data from the given
			%						sample where N is the number of cells in the
			%						returned data and M is the number of
			%						channels requested. 
			%
			%		gateLogicals	<cell> S x G cell array of N x 1 logical vector 
			%						indicating which points of sample ID S are in the 
			%						given gates (G) and/or bins. 
			
			% Check and extract inputs
			[sliceData, sliceChannels, sliceDataType, sliceGates] = zCheckInputs_slice(self);
			
			% Slice out data
			dataMatrix = [];
			gateLogicals = cell(size(sliceGates));
			for s = 1:numel(sampleIDs)
				
				if isempty(sliceGates{s, 1}) % Some tcData will be empty 
					continue
				end
				
				% Extract data
				sID = sampleIDs(s);
				dataS = zeros(numel(sliceGates{s, 1}), numel(sliceChannels));
				for ch = 1:numel(sliceChannels)
					if ismember(sliceChannels, Gating.SCATTER_CHANNELS)
						sdt = 'raw';
					else
						sdt = sliceDataType;
					end
					dataS(:, ch) = sliceData(sID).(sliceChannels{ch}).(sdt)(sliceGates{s, ch});
					
					gateLogicals{s, ch} = false(sliceData(sID).nObs, 1);
					gateLogicals{s, ch}(sliceGates{s, ch}) = true;
				end
				dataMatrix = [dataMatrix; dataS];
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
				chanFieldnames = {'channel', 'Channel', 'channels', 'Channels'};
				doChan = isfield(sliceParams, chanFieldnames);
				if any(doChan)
					chanF = chanFieldnames{doChan};
					sliceChannels = sliceParams.(chanF);
					
					validateattributes(sliceChannels, {'char', 'cell'}, {}, mfilename, ['sliceParams.', chanF]);
					if ischar(sliceChannels), sliceChannels = {sliceChannels}; end % Force cell
					sliceChannels = reshape(sliceChannels, 1, []); % Force row vector
					
					badChannels = setdiff(sliceChannels, [self.channels, Gating.SCATTER_CHANNELS, 'Time']);
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
				sampleIDs = reshape(unique(round(sampleIDs)), 1, []); 
				assert(all(sampleIDs <= numel(sliceData)), ...
					'At least one sampleID is too large')
				
				% Slices data of the given dataType
				if isfield(sliceParams, 'dataType')
					validatestring(sliceParams.dataType, self.dataTypes, mfilename, 'sliceParams.dataType');
					sliceDataType = sliceParams.dataType;
				else
					sliceDataType = 'raw'; % Default is raw data
				end
				
				% Applies the given gate to the sliced data
				gateFieldnames = {'gate', 'Gate', 'gates', 'Gates'};
				doGate = isfield(sliceParams, gateFieldnames);
				if any(doGate)
					gateF = gateFieldnames{doGate};
					gates = sliceParams.(gateF);
					
					validateattributes(gates, {'char', 'cell'}, {}, mfilename, ['sliceParams.', gateF]);
					if ischar(gates), gates = {gates}; end
					assert(numel(gates) == 1 | numel(gates) == numel(sliceChannels), ...
							'Number of gates must be 1 or match number of channels');
					sliceGates = cell(numel(sampleIDs), numel(gates));
					
					for gi = 1:numel(gates)
						currGate = gates{gi};
						
						% Flip gates with 'x' in front that don't already exist as 'xABC' gates
						if strcmp(currGate(1), 'x') && ~ismember(currGate, self.gateNames)
							currGate = currGate(2:end);
							flipGate = true;
						else
							flipGate = false;
						end
						
						validatestring(currGate, self.gateNames, mfilename, ['sliceParams.' gateF, '{', num2str(gi), '}']);
						
						for si = 1:numel(sampleIDs)
							sampleID = sampleIDs(si);
							if ~isempty(sliceData(sampleID).(sliceChannels{1})) % Some tcData will be empty 
								sliceGate = sliceData(sampleID).gates.(currGate);
								if flipGate, sliceGate = ~sliceGate; end
		% 						sliceGates{si} = sliceGate; % old code for logical indexing
								sliceGates{si, gi} = find(sliceGate);
							end
						end
					end
				else
					sliceGates = cell(numel(sampleIDs), 1);
					for si = 1:numel(sampleIDs)
						sampleID = sampleIDs(si);
						if ~isempty(sliceData(sampleID).(sliceChannels{1})) % Some tcData will be empty 
							numSliceCells = numel(sliceData(sampleID).(sliceChannels{1}).raw);
	% 						sliceGates{si} = true(numSliceCells, 1); % old code for logical indexing
							sliceGates{si} = (1:numSliceCells)'; % Default is all cells
						end
					end
				end
				% If just one gate given, expand to cover all channels 
				if ((size(sliceGates, 2) == 1) && (numel(sliceChannels) > 1))
					sliceGates = repmat(sliceGates, 1, numel(sliceChannels));
				end
				
				% Find minimal number of points to take from each sample
				pointsPerGate = cellfun(@numel, sliceGates);
				if (isfield(sliceParams, 'numPoints') && ~isempty(sliceParams.numPoints))
					% Determine if a minimal number of points have been set
					numPoints = sliceParams.numPoints;
				else
					numPoints = inf;
				end
				if (isfield(sliceParams, 'equalize') && sliceParams.equalize)
					% Use the same number of points per sample
					numPoints = min(numPoints, min(min(pointsPerGate)));
					
					for sgi = 1:numel(sliceGates) % Do linearly because multiple dimensions
						if ~isempty(sliceGates{sgi}) % Some tcData will be empty 
							ss = FlowAnalysis.subSample(numel(sliceGates{sgi}), numPoints);
							sliceGates{sgi} = sliceGates{sgi}(ss);
						end
					end
				else
					% Use the same number of points per gate per sample
					numPoints = min(numPoints, min(pointsPerGate, [], 2));
					
					for si = 1:size(sliceGates, 1)
						for gi = 1:size(sliceGates, 2)
							if ~isempty(sliceGates{si, gi}) % Some tcData will be empty 
								ss = FlowAnalysis.subSample(numel(sliceGates{si, gi}), numPoints(si));
								sliceGates{si, gi} = sliceGates{si, gi}(ss);
							end
						end
					end
				end
				
				
				% Slices from the given bins. Controls are not binned, so if they 
				% are the sliceData, we just skip the rest of the inputs checking
				if (isfield(sliceParams, 'bins') && ~isempty(self.bins) && ~sliceControls)
					
					validateattributes(sliceParams.bins, {'numeric', 'char'}, {}, mfilename, 'sliceParams.bins');
					if ischar(sliceParams.bins)
						if strcmpi(sliceParams.bins, 'all')
							sliceParams.bins = (1:self.numBins)';
						else
							sliceParams.bins = []; 
						end
					end
				else
					sliceParams.bins = [];
				end
				
				% Only bin if we have some bins to bin!
				if ~isempty(sliceParams.bins)
					assert(size(sliceParams.bins, 2) == 1 || ...
						   size(sliceParams.bins, 2) == numel(fieldnames(self.binInputs)), ...
						   'Bin IDs formatted incorrectly')
					assert(size(sliceParams.bins, 1) <= self.numBins, ...
						   'Too many bins requested')

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
					for si = 1:numel(sampleIDs)
						sampleID = sampleIDs(si);
						cellsInBins = [];
						for b = binIdxs % Sequentially looks at each bin ID
							cellsInBins = [cellsInBins, self.bins{sampleID}{b}];
						end
						for gi = 1:size(sliceGates, 2) % Handles gates for each channel
							if ~isempty(sliceGates{si, gi}) % Some tcData will be empty 
								sliceGates{si, gi} = intersect(sliceGates{si, gi}, cellsInBins);
							end
						end
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
					'Incorrect scalar length passed');
				
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
				
		
		function self = threshGate(self, channels, mode, thresh, options, sampleIDs)
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
			%	self = self.threshGate(channels, mode, thresh, options)
			%
			%	Inputs
			%		channels	<char, cell> A cell list of channels to threshold on. 
			%					(A single string for one channel is also accepted)
			%
			%		mode		(optional) <char> {'or', 'and'} - determines how 
			%					thresh gates are crossed (if multiple channels given)
			%
			%		thresh		(optional) Either a specific threshold value to use 
			%					or one of the keywords: {'manual', 'auto'}. In
			%					any case, the threshold value and any generated
			%					plots are saved in the 'Gating' folder. 
			%						Numeric: The thresh value given should be a
			%						non-transformed number. A single value can be 
			%						given to threshold in each channel, or individual 
			%						values for each channel. 
			%						'Manual': The user will draw a line on a plot, 
			%						the intercept of which determines the threshold
			%						 --> Default if no option given
			%						'Auto': The thresh is automatically determined
			%						from the 99.9th %ile of the wtData. 
			%
			%		options		(optional) <char, cell> A string or cell list of 
			%					strings with additional options: 
			%						'recompute': Overwrites the automatic loading 
			%						of existing threshold values from saved files. 
			%						<'dataType'> Uses the specified dataType for 
			%						defining the thresholds. 
			%
			%		sampleIDs	(optional) <numeric> A subset of sampleIDs to use 
			%					for generating the thresholds. By default, all the
			%					samples are used.
			
			zCheckInputs_threshGate(self);
			
			% Setup new directory for gates
			gateDir = [self.folder, 'Gating', filesep];
			if ~exist(gateDir, 'file')
				mkdir(gateDir)
			end
			gateSaveName = [self.date, '_', self.name];
			threshFname = [gateDir, 'Thresh-Vals_', gateSaveName, '.mat'];
			
			% Load thresholds from file or from the object's fields if possible
			if ~ismember('recompute', options)
				if ~isempty(self.threshGateVals)
					threshVals = self.threshGateVals;
				elseif exist(threshFname, 'file')
					load(threshFname, 'threshVals');
				end
			end
			
			% Check which dataType to use for thresholding
			if any(ismember(self.dataTypes, options))
				dti = find(ismember(self.dataTypes, options), 1);
				dtype = self.dataTypes{dti};
			else
				dtype = 'raw';
			end
			
			% If there were no thresh values to load, create them
			if ~exist('threshVals', 'var') 
				threshVals = zeros(1, numel(channels));
				for ci = 1:numel(channels)
					
					chan = channels{ci};
					
					if (ischar(thresh) && strcmpi(thresh, 'auto'))
						% Thresh value is found using the untransfected (wt) cells
						% --> Set at 99.9th %ile of fluorescence in a channel
						threshVals(ci) = prctile(self.controlData(end).(chan).(dtype), 99.9);
					elseif isnumeric(thresh)
						% Thresh value given as an option
						threshVals(ci) = thresh(ci); % Ignore extra entries
					else
						% Thresh value is found by marking a line on a graph
						combData = [];
						for i = sampleIDs
							combData = [combData; self.sampleData(i).(chan).(dtype)(1:self.numSamples:end)];
						end
						
						pos = find(combData > 0);
						ss = FlowAnalysis.subSample(numel(pos), 1e4);
						combDataLog = log10(combData(pos(ss)));
						
						figThresh = figure();
						ax = gca(); hold(ax, 'on')
						histogram(ax, combDataLog)
						title('Draw a line to set an x-axis threshold', 'fontsize', 16)
						ylabel('Count', 'fontsize', 14)
						xlabel(['log_{10} ', strrep(chan, '_', '-')], 'fontsize', 14)
						
						h = imline();
						position = wait(h);
						self.gatePolygons.(['TH', self.SHORT_COLORS.(chan)]) = position;
						
						saveas(figThresh, [gateDir, 'Gate-Thresh_', ...
								strrep(chan, '_', '-'), '_', gateSaveName, '.fig']);
						
						% Fit the points to get a line, then use that to find the
						% x-intercept (threshold value)
						fit = polyfit(position(:, 1), position(:, 2), 1);
						threshVals(ci) = 10.^(-fit(2) / fit(1));
					end
				end
			end
			
			% Save thresh values
			save(threshFname, 'threshVals');
			self.threshGateVals = threshVals; 
			
			% Apply thresh gate
			thrGateNames = cell(1, numel(channels));
			for ci = 1:numel(channels)
				
				chan = channels{ci};
				
				% Find cells passing threshold and record gate
				thrGateNames{ci} = ['TH', self.SHORT_COLORS.(chan)];
				for i = 1:self.numSamples
					passThresh = (self.sampleData(i).(chan).(dtype) >= threshVals(ci));
					self.sampleData(i).gates.(thrGateNames{ci}) = passThresh;
				end
				
				% Do the same for controls
				for i = 1:self.numControls
					if ~isempty(self.controlData(i).(chan)) % Some tcData will be empty 
						passThresh = (self.controlData(i).(chan).(dtype) >= threshVals(ci));
						self.controlData(i).gates.(thrGateNames{ci}) = passThresh;
					end
				end
			end
			
			self.addGates(thrGateNames);
			if exist('mode', 'var') && ~isempty(mode) && (numel(channels) > 1)
				self.crossGates(thrGateNames, mode);
			end
			
			fprintf(1, 'Finished thresholding\n');
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_threshGate(self)
				validateattributes(channels, {'char', 'cell'}, {}, mfilename, 'gates', 1);
				if ischar(channels), channels = {channels}; end % For simplicity
				badChannels = setdiff(channels, self.channels);
				assert(isempty(badChannels), 'Channel not valid: %s\n', badChannels{:});

				if exist('mode', 'var') && ~isempty(mode) && (numel(channels) > 1)
					mode = lower(mode);
					validatestring(mode, {'and', 'or'}, mfilename, 'mode', 2);
				end
				
				if exist('thresh', 'var')
					validateattributes(thresh, {'char', 'numeric'}, {}, mfilename, 'thresh', 3);
					if isnumeric(thresh)
						if (numel(thresh) == 1)
							thresh = repmat(thresh, size(channels));
						end
						assert(numel(thresh) == numel(channels), ...
							'# threshold values given different than # of channels given');
					else
						thresh = lower(thresh);
					end
				else
					thresh = false;
				end
				
				if exist('options', 'var')
					validateattributes(options, {'cell', 'char'}, {}, mfilename, 'options', 4);
					if ischar(options), options = {options}; end % For simplicity
					options = cellfun(@lower, options, 'uniformoutput', false);
				else
					options = {};
				end
				
				if exist('sampleIDs', 'var')
					validateattributes(sampleIDs, {'numeric'}, {'positive'}, mfilename, 'sampleIDs', 5);
					assert(max(sampleIDs) <= self.numSamples, 'Sample ID(s) too large');
				else
					sampleIDs = 1:self.numSamples;
				end
			end
		end
				
		
		function self = crossGates(self, gates, mode)
			% Crosses the given gates using the given crossing mode
			% A new gate is created with name in the following form: 
			%		'gate1_gate2_gate3_[...]_gateN'
			%
			%	self = self.crossGates(gates, mode)
			%
			%	Inputs
			%		gates	<cell> A cell list of gate names to cross. Must be >= 2 gates
			%		mode	<char> {'or', 'and'} - determines how gates are crossed
			
			% Check inputs
			validateattributes(gates, {'cell'}, {'vector'}, mfilename, 'gates', 1);
			validatestring(mode, {'and', 'or'}, mfilename, 'mode', 2);
			newGateName = [gates{1}, sprintf('_%s', gates{2:end})];
			
			flipGates = false(1, numel(gates));
			for gi = 1:numel(gates)
				% Flip gates with 'x' in front that don't already exist as 'xABC' gates
				if strcmp(gates{gi}(1), 'x') && ~ismember(gates{gi}, self.gateNames)
					gates{gi} = gates{gi}(2:end);
					flipGates(gi) = true;
				end
			end
			
			badGates = setdiff(gates, self.gateNames);
			assert(isempty(badGates), 'Gate does not exist: %s\n', badGates{:});
			
			% Do gate crossing
			switch mode
				case 'or'
					for i = 1:self.numSamples
						self.sampleData(i) = crossOR(self.sampleData(i), newGateName, flipGates);
					end
					for i = 1:self.numControls
						if ~isempty(self.controlData(i).(self.channels{1})) % Some tcData will be empty 
							self.controlData(i) = crossOR(self.controlData(i), newGateName, flipGates);
						end
					end
				case 'and'
					for i = 1:self.numSamples
						self.sampleData(i) = crossAND(self.sampleData(i), newGateName, flipGates);
					end
					for i = 1:self.numControls
						if ~isempty(self.controlData(i).(self.channels{1})) % Some tcData will be empty 
							self.controlData(i) = crossAND(self.controlData(i), newGateName, flipGates);
						end
					end
			end
			
			self.addGates(newGateName);
			fprintf(1, 'Finished crossing gates\n')
			
			
			% --- Helper functions --- %
			
			
			function data = crossOR(data, newGateName, flipGates)
				% Crosses w/ OR logic
				
				crossedGates = false(size(data.gates.(gates{1})));
				for g = 1:numel(gates)
					currGate = data.gates.(gates{g});
					if flipGates(g), currGate = ~currGate; end
					crossedGates = (crossedGates | currGate);
				end
				data.gates.(newGateName) = crossedGates;
			end
			
			
			function data = crossAND(data, newGateName, flipGates)
				% Crosses w/ AND logic 
				
				crossedGates = true(size(data.gates.(gates{1})));
				for g = 1:numel(gates)
					currGate = data.gates.(gates{g});
					if flipGates(g), currGate = ~currGate; end
					crossedGates = (crossedGates & currGate);
				end
				data.gates.(newGateName) = crossedGates;
			end
		end
		
		
		function [gatePcts] = printGatePcts(self)
			% Prints the percent of cells passing each defined gate into a table
			%
			%	gatePcts = self.printGatePcts();
			%
			%	Outputs
			%		gatePcts	<table> Gate percentages in a table format.
			%					Table headers are gate names. 
			
			% Some gates may not be applied to the sample data, so only look at
			% names of gates literally assigned under sampleData
			sampleGateNames = fieldnames(self.sampleData(1).gates);
			numGates = numel(sampleGateNames);
			gatePcts = zeros(self.numSamples, numGates);
			
			for g = 1:numGates
				gate = sampleGateNames{g};
				
				for si = 1:self.numSamples
					% Use mean as simple way to compute fraction passing gate
					gatePcts(si, g) = mean(self.sampleData(si).gates.(gate)) * 100;
				end
			end
			
			% Convert to readable table
			sampleGateNames = fieldnames(self.sampleData(1).gates);
			gatePcts = array2table([self.numCells, gatePcts]);
			gatePcts.Properties.VariableNames = [{'Num_cells'}; sampleGateNames];
			
			% Append to Sample Map
			gatePcts = [self.sampleMap, gatePcts];
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
				assert(logical(exist(newSampleMap, 'file')), 'New Sample Map file not found');
				newSampleMap = readtable(newSampleMap, 'delimiter', '\t');
			elseif iscell(newSampleMap)
				assert(self.sampleMap.width == size(newSampleMap, 2), ...
					'New Sample Map has incorrect number of variables (%d) compared to the data (%d)', ...
					size(newSampleMap, 2), self.sampleMap.width)
				varNames = self.sampleMap.Properties.VariableNames;
				newSampleMap = cell2table(newSampleMap);
				newSampleMap.Properties.VariableNames = varNames;
			end
			
			% Adjust table size for combined replicates
			if ismember('Replicate', newSampleMap.Properties.VariableNames)
				% Wean sampleMap to only be have one replicate
				newSampleMap = newSampleMap(newSampleMap.Replicate == 1, :);
			end
			
			% Check table size
			assert(height(newSampleMap) == self.numSamples, ...
				'New Sample Map has incorrect number of samples (%d) compared to the data (%d)', ...
				height(newSampleMap), self.numSamples);
			
			% Re-assign variable
			self.sampleMap = newSampleMap;
			
			fprintf(1, 'Finished adding new sample map\n');
		end
		
		
		function self = editLogicleParams(self, newParams)
			% Allows the user to change the logicle conversion parameters used
			% by the object by passing a new set of parameters. This method
			% ensures that the new parameters are valid.
			
			validateattributes(newParams, {'struct'}, {}, mfilename, 'newParams', 1);
			
			% Look for any missing fields - keep existing parameters for those
			% that are not given
			if ~isfield(newParams, 'T'), newParams.T = self.logicleParams.T; end
			if ~isfield(newParams, 'M'), newParams.M = self.logicleParams.M; end
			if ~isfield(newParams, 'r'), newParams.r = self.logicleParams.r; end
			
			% Check the values are valid
			validateattributes(newParams.T, {'numeric'}, {'scalar', 'positive'}, mfilename, 'newParams.T');
			validateattributes(newParams.M, {'numeric'}, {'scalar', 'positive'}, mfilename, 'newParams.M');
			validateattributes(newParams.r, {'numeric'}, {'scalar', 'negative'}, mfilename, 'newParams.r');
			
			self.logicleParams = newParams;
		end
		
		
		function rot = rotate(self, channels, theta)
			% Rotates the data by the given angle (theta, degrees), returning a new
			% FlowData object with rotated data, but the same coordinate system
			%
			%	Inputs
			%		channels	<char, cell> Two channels to rotate on
			%
			%		theta		<numeric> The angle to rotate by. Negative = clockwise
			%
			%	Outputs
			%		rot			<FlowData> A handle to the new object w/ rotated data
			%
			%	
			%	Additional information: 
			%
			%		The rotation matrix:
			%
			%			[cos(theta), -sin(theta)
			%			 sin(theta),  cos(theta)]
			%
			%		is multiplied against a two-row matrix where the X-values
			%		are in the first row and the Y-values are in the second row.
			%		The points are then rotated theta degress counter-clockwise
			%		about the origin. 
			%
			%		<from: https://en.wikipedia.org/wiki/Rotation_matrix >
			
			% Check angle input
			zCheckInputs_rotate(self);
			
			% Define rotation matrix
			thetaRads = theta * pi / 180; 
			rotMatrix = [cos(thetaRads), -sin(thetaRads)
						 sin(thetaRads),  cos(thetaRads)];
			
			% Create output data structure
			rot = self.copy();
			
			% Rotate data
			for di = 1:numel(self.dataTypes)
				
				dt = self.dataTypes{di};
				
				for si = 1:self.numSamples
					
					% Extract data
					sd = [rot.sampleData(si).(channels{1}).(dt), ...
						  rot.sampleData(si).(channels{2}).(dt)]';
					
					% Apply rotation
					sd = rotMatrix * sd;
					
					% Save data to new data structure
					rot.sampleData(si).(channels{1}).(dt) = sd(1, :)';
					rot.sampleData(si).(channels{2}).(dt) = sd(2, :)';
				end
				
% 				for ci = 1:self.numControls
% 					
% 					% Extract data
% 					cd = [rot.controlData(ci).(channels{1}).(dt), ...
% 						  rot.controlData(ci).(channels{2}).(dt)]';
% 					
% 					% Apply rotation
% 					cd = rotMatrix * cd;
% 					
% 					% Save data to new data structure
% 					rot.controlData(ci).(channels{1}).(dt) = cd(1, :)';
% 					rot.controlData(ci).(channels{2}).(dt) = cd(2, :)';
% 				end
			end
			
			
			% --- Helper Functions --- %
			
			
			function zCheckInputs_rotate(self)
				
				validateattributes(channels, {'cell', 'char'}, {}, mfilename, 'channels', 1);
				validateattributes(theta, {'numeric'}, {'scalar', 'real'}, mfilename, 'theta', 2);
				
				badChannels = setdiff(channels, self.channels);
				assert(isempty(badChannels), 'Channel not recognized: %s', badChannels{:});
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
			for gi = 1:numel(gates)
				if ~any(strcmp(gates{gi}, self.gateNames))
					self.gateNames = [self.gateNames, gates(gi)];
					added = [added, gates(gi)];
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