function [PFit, minVal, fitPDF, fitPDF_log, ubID] = fitPdfTX3(xdata, options)
	% Fits the transfection distribution to a sample x
	%	
	%	[PFit, minVal, fitPDF, fitPDF_log, ubID] = ...
	%		fitPdfTX(xdata, gamOnly, numIter, doMEF, logicleParams)
	%	
	%	PFit array order: ['k', 'th', 'mu', 'sigma']
	%
	%	options fields: 'edgesLin', 'edgesLog', 'targetsLin', 'P0', 
	%					'numRep', 'numIter', 'scale', 'logicleParams'
	%
	%	Note: density is logicle-transformed to ensure that the whole distribution 
	%	is weighted equally. 
	%
	%		Alt: apply a weight vector which is exponentially-decreasing
	
	% Check inputs
	zCheckInputs_fitPdfTX3();
	
	% Find density of data in linear and logicle space
	densityLin = histcounts(xdata, options.edgesLin, 'Normalization', 'Probability');
	densityLog = histcounts(Transforms.lin2logicle(xdata, true, options.logicleParams), ...
			options.edgesLog, 'Normalization', 'Probability');
	centersLin = options.edgesLin(2:end) - diff(options.edgesLin) / 2;
	centersLog = options.edgesLog(2:end) - diff(options.edgesLog) / 2;
	
	% Remove Inf and NaN values
	validLin = (densityLin > 0 & ~isinf(densityLin) & ~isnan(densityLin));
	validLog = (densityLog > 0 & ~isinf(densityLog) & ~isnan(densityLog));
	
	% Find location of linear bins within logicle bins
	[bIDs, ubID] = findBins(options.edgesLog, centersLin(validLin));
	
	minVals = inf(1, options.numRep);
	P = repmat(options.P0, [1, options.numRep]);
	for ri = 1:options.numRep
		
		dP = ones(size(options.P0));
		
		for i = 1:options.numIter
			% dP is a change vector applied to the parameters
			%	It initizializes as 1, then changes according to whether 
			%	it makes the fit better or worse. Better = repeat change again
			%	Worse = try something else. 
			P0_i = P(:, ri) .* dP;
			fval = minfunc(P0_i, bIDs, ubID);
			if (fval < minVals(ri))
				minVals(ri) = fval;
				P(:, ri) = P0_i;
			else
				dP(:) = 1;
				pi = randi(numel(options.P0));
				dP(pi) = 10.^(randn(1) / 5);
% 				dP =  10.^(randn(size(dP)) / 5);
			end
		end
	end
% 	minVals, P
	[minVal, minIdx] = min(minVals);
	PFit = P(:, minIdx);
	
	% Find parameters
% 	options = optimset('Display', 'off');
% 	[P, fval] = fminsearch(@minfunc, P0, options);
	
	% Generate final PDF
	[bIDs, ubID] = findBins(options.edgesLog, options.targetsLin);
	fitPDF = pdfTX(options.targetsLin, PFit(1), PFit(2), PFit(3), PFit(4));
	fitPDF_log = makeLogPDF(fitPDF, bIDs);
	
	fprintf(1, 'Fval = %.2f\n', minVal);
	
	
	% --- Helper Functions --- %
	
	
	function [bIDs, ubID] = findBins(logEdges, linCenters)
		
		[~, ~, bID] = histcounts(Transforms.lin2logicle( ...
				linCenters, true, options.logicleParams), logEdges);
		ubID = unique(bID(bID ~= 0));
		bIDs = cell(1, numel(ubID));
		for ubi = 1:numel(ubID)
			bIDs{ubi} = (bID == ubID(ubi));
		end
	end
	
	
	function fval = minfunc(p, bIDs, ubID)
		testPDF = pdfTX(centersLin, p(1), p(2), p(3), p(4));
% 		xdiff = numel(xvals_lin) - numel(testPDF); % Needed when min(xvals) > 0 
		
		switch options.mode
			case {'log', 'logicle'}
				logTestPDF = makeLogPDF(testPDF, bIDs);
				fval = sum(abs(log10(logTestPDF) - log10(densityLog(ubID))).^2);
			case {'lin', 'linear'}
				fval = sum(abs(log10(testPDF(validLin)) - log10(densityLin(validLin))).^2);
			% Should only have to find valid points in the lin data since the
			% log bins are determined based on valid lin data already
		end
	end


	function logPDF = makeLogPDF(linPDF, bIDs)
		logPDF = zeros(size(bIDs));
		for bi = 1:numel(bIDs)
			logPDF(bi) = sum(linPDF(bIDs{bi}));
		end
	end


	function zCheckInputs_fitPdfTX3()
		
		% Logicle parameters
		if isfield(options, 'scale')
			validateattributes(options.scale, {'numeric'}, {'scalar', 'positive'}, mfilename, 'scale');
		else
			options.scale = 1; 
		end
		if isfield(options, 'logicleParams')
			validateattributes(options.logicleParams, {'struct'}, {}, mfilename, 'logicleParams');
			if isfield(options.logicleParams, 'MEF')
				% Overwrite scale if it is supplied in logicleParams
				options.scale = options.logicleParams.MEF; 
			end
		else
			options.logicleParams = struct('MEF', options.scale);
		end
		options.logicleParams = Transforms.checkLogicleParams(true, options.logicleParams);
		
		% X values and bins
		validateattributes(xdata, {'numeric'}, {'vector'}, mfilename, 'xdata', 1);
		if isfield(options, 'edgesLin')
			validateattributes(options.edgesLin, {'numeric'}, {'vector'}, mfilename, 'edgesLin');
		else
			% Standard AFU range adjusted for logicle scale
			options.edgesLin = -1e2:1e2:1e5 * options.scale;
		end
		if isfield(options, 'edgesLog')
			validateattributes(options.edgesLog, {'numeric'}, {'vector'}, mfilename, 'edgesLog');
		else
			% Standard logicle-transformed space (indpendent of scale)
			options.edgesLog = 0:0.1:4.5;
		end
		if isfield(options, 'targetsLin')
			validateattributes(options.targetsLin, {'numeric'}, {'vector'}, mfilename, 'targetsLin');
		else
			options.targetsLin = options.centersLin;
		end
		
		% Parameters
		if isfield(options, 'P0')
			validateattributes(options.P0, {'numeric'}, {'vector', 'numel', 4}, mfilename, 'P0');
		else
			% Initial values for [k, th, mu, sigma]
			options.P0 = [0.1; 3e3 * options.scale; 0; 30 * options.scale];
		end
		
		% Fit Replicates/Iterations
		if isfield(options, 'numRep')
			validateattributes(options.numRep, {'numeric'}, {'scalar', 'positive'}, mfilename, 'numRep');
		else
			options.numRep = 5;
		end
		if isfield(options, 'numIter')
			validateattributes(options.numIter, {'numeric'}, {'scalar', 'positive'}, mfilename, 'numIter');
		else
			options.numIter = 1e3;
		end
		
		% Fit mode
		if isfield(options, 'mode')
			validatestring(options.mode, {'lin', 'linear', 'log', 'logicle'}, mfilename, 'mode');
		else
			options.mode = 'log';
		end
		
	end
	
end


