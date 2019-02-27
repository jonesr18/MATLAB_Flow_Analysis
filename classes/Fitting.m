classdef Fitting < handle
	
	methods (Static)
		
		function [newPdf, binIdxs, newXvals] = adjustPdfBins(oldPdf, oldXvals, newEdges)
			% Adjusts a PDF by summing the density within newly defined bins
			%
			%	[newPdf, newXvals] = adjustPdfBins(xvals, pdf, newXedges);
			%
			%	Inputs
			%
			%		oldPdf		<numeric> The probability density function (PDF)
			%					to be adjusted
			%
			%		oldXvals	<numeric> The x-values of the points in the
			%					starting PDF
			%
			%		newEdges	<numeric> The x-values of the edges of the bins
			%					used to create the new PDF. 
			%
			%	Outputs
			%
			%		newPdf		<numeric> The newly generated PDF
			%
			%		binIdxs		<numeric> The indexes corresponding to valid bins
			%
			%		newXvals	<numeric> The new xvals for the new PDF. Computed 
			%					automatically by averaging each pair of edges in 
			%					newXedges.

			
			
			% Find the location of each starting xval in the newEdges
			%	binAssignments has the assigned bin ID for each oldXval
			[~, ~, binAssignments] = histcounts(oldXvals, newEdges);
			numNewBins = numel(newEdges) - 1;
			subBinIdxs = cell(1, numNewBins);
			for bi = 1:numNewBins
				subBinIdxs{bi} = (binAssignments == bi);
			end
			
			% Find the new bins which actually have points
			binIdxs = unique(binAssignments(binAssignments ~= 0));
			
			% Sum all oldPdf values within new bins
			newPdf = zeros(size(subBinIdxs));
			for sbi = 1:numel(subBinIdxs)
				newPdf(sbi) = sum(oldPdf(subBinIdxs{sbi}));
			end
			
			% Find new x values by averaging new edges (unless pre-supplied)
			newXvals = zeros(1, numNewBins);
			for ei = 1:numNewBins
				newXvals(ei) = mean(newEdges([ei, ei + 1]));
			end
		end
		
		
		function [centers, edges, centers2, c2_idx] = createSpace(xmin, xmax, dx)
			% Generates centers, edges (for binning), the edges that result from
			% convolution of two centers-based pdf vectors, and the indexes
			% corresponding with centers in centers2. 

			% TODO: should xmin shift right when convoluting xmin > 0?

			assert(xmin < xmax, 'xmin must be less than xmax');
			assert(dx < (xmax - xmin), 'dx too large');

			edges = (xmin - dx / 2) : dx : (xmax + dx / 2);
			centers = xmin : dx : xmax;

			if sign(xmin) ~= sign(xmax)
				centers2 = ((2 * xmin) : dx : (2 * xmax)) - sign(xmin) * mod(xmin, dx);
			else % No changes made yet
				centers2 = ((2 * xmin) : dx : (2 * xmax)) - sign(xmin) * mod(xmin, dx);
			end
			
			% The centers need to be rounded to avoid floating-point errors
			%	-1 * (...) since rounding decimals are opposite log10 vals
			roundDec = -1 * (floor(log10(dx)) - 1);		
			c2_idx = ismember(round(centers2, roundDec), round(centers, roundDec));
		end
		
		
		function xpdf = fillPdfInfs(xpdf, dx)
			% Replaces inf and nan values in a pdf with the value that should be there, 
			% assuming the entire pdf * dx sums to 1

			x0 = (isinf(xpdf) | isnan(xpdf));
			p0 = (1 - sum(xpdf(~x0)) .* dx) ./ dx;	% Adjust so sum * dx = 1
			xpdf(x0) = p0 ./ sum(x0);

		end
		
		
		function [gapdf, gpdf, npdf] = pdfGammaNormal(x, k, theta, mu, sigma, C)

			assert(numel(unique(round(diff(x)))) == 1, 'x must be evenly spaced for convolution')
			% Note, the above check needs 'round' due to floating point error
			if exist('C', 'var')
				assert(size(C, 1) == size(C, 2), 'C must be square!'); 
			else
				C = 1; % For co-variance - but not used currently
			end

			gpdf = gampdf(x, k, theta);
		% 	gpdf = gampdf_ND(x, k, theta, C);
			npdf = normpdf(x, mu, sigma);	% No covariance in autofluorescence

			% Find x step value and generate x vals for conv vector
			dx = mean(diff(x));
			[~, ~, ~, x2_idx] = Fitting.createSpace(min(x), max(x), dx);

			% Gamma pdf = inf at x = 0, so take that out and replace it w/ 
			% remaining probability from the summation of the rest
			gpdf = Fitting.fillPdfInfs(gpdf, dx);

			% Convolve and extract the conserved x values
			% --> Convolution basically extends the x vector 2*x in both directions from
			% 0, but since there is no x vector given, we have to adjust x manually.
			gapdf = fconv(gpdf, npdf);
			gapdf = gapdf(x2_idx) .* dx;

			% Normalize
			gapdf = gapdf ./ sum(gapdf);
		end
		
		
		function [gapdf, ppdf, gpdf, npdf] = pdfPoissGammaNormal(x, lambda, k, theta, mu, sigma)

			assert(numel(unique(diff(x))) == 1, 'x must be evenly spaced for convolution')

			xpoiss = 0:100;
			ppdf = poisspdf(xpoiss, lambda); 
			gpdf = zeros(size(x)); 
			for pi = 1:numel(ppdf)
				gpdf = gpdf + ppdf(pi) * gampdf(x, k * xpoiss(pi), theta);
			end
			npdf = normpdf(x, mu, sigma);

			% Find x step value and generate x vals for conv vector
			dx = mean(diff(x));
			[~, ~, ~, x2_idx] = Fitting.createSpace(min(x), max(x), dx);

			% Gamma pdf = inf at x = 0, so take that out and replace it w/ 
			% remaining probability from the summation of the rest
			gpdf = Fitting.fillPdfInfs(gpdf, dx);

			% Convolve and extract the conserved x values
			% --> Convolution basically extends the x vector 2*x in both directions from
			% 0, but since there is no x vector given, we have to adjust x manually.
			gapdf = fconv(gpdf, npdf);
			gapdf = gapdf(x2_idx) .* dx;

			% Normalize
			gapdf = gapdf ./ sum(gapdf);
		end
		
		
		function [rappdf, rpdf, npdf, ppdf] = pdfPoissRadiiNormal(x, lambda, mu1, sig1, mu2, sig2)
			
			% TODO - this isn't working
			
			assert(numel(unique(diff(x))) == 1, 'x must be evenly spaced for convolution')
			
			% Find x step value and generate x vals for conv vector
			dx = mean(diff(x));
			[~, ~, ~, x2_idx] = Fitting.createSpace(min(x), max(x), dx);
			
			% Multiply Poisson dist by radii dist by convolving radiiPdf by
			% itself over and over, weighting by the probability of each point
			% in the poisson distribution. 
			xpoiss = 0 : (lambda * 100);
			ppdf = poisspdf(xpoiss, lambda); 
			rpdf = zeros(size(x));
			rpdf(x == 0) = ppdf(1);
			for i = 2:numel(ppdf)
				% This func is inf at 0, so fix that before doing anything
				rpdfTemp = Fitting.fillPdfInfs(radiiPdf(x, mu1, sig1), dx);
				rpdfNew = rpdfTemp;
				for j = 2:xpoiss(i) % Do only for 2+ complexes/cell
					rpdfTemp = fconv(rpdfTemp, rpdfNew);
					rpdfTemp = rpdfTemp(x2_idx) .* dx;
				end
				rpdf = rpdf + ppdf(i) .* rpdfTemp ./ sum(rpdfTemp); 
			end
			
			% Convolve and extract the conserved x values
			% --> Convolution basically extends the x vector 2*x in both directions from
			% 0, but since there is no x vector given, we have to adjust x manually.
			npdf = normpdf(x, mu2, sig2);
			rappdf = fconv(rpdf, npdf ./ sum(npdf));
			rappdf = rappdf(x2_idx) .* dx;

			% Normalize
			if any(isinf(rappdf) | isnan(rappdf))
				warning('PDF contains infs or nans!')
			end
			rappdf = rappdf ./ sum(rappdf);
			
			
			% --- Helper Function --- %
			
			
			function rpdf = radiiPdf(x, mu, sig)
				rpdf = 1 / 3 ./ x .* (3 / 4 / pi .* x).^(1/3) .* ...
						(2 * pi * sig^2)^(-1/2) .* ...
						exp(-1 / (2 * sig^2) * ((3 / 4 / pi * x).^(1/3) - mu).^2); 
			end
			
			
		end
		
		
		function protPerCell = radiiRnd(numCells, lambda, mu1, sig1, sf, mu2, sig2, numPoly)
			
			if ~exist('numPoly', 'var')
				numPoly = 1;
			end
			
			% Poisson for complexes per cell
			% See link below - sum of poissons = poisson w/ summed lambda
			% https://en.wikipedia.org/wiki/List_of_convolutions_of_probability_distributions
			complexesPerCell = poissrnd(lambda, numCells, 1);
			
			% Find amount of DNA (in terms of radii) per complex, then convert
			% to volume by relating radius to volume of a sphere
			dnaPerCell = zeros(numCells, numPoly);
			for i = 1:numCells
				complexAssignment = randi(numPoly, 1, complexesPerCell(i));
				for j = 1:numPoly
					complexRadii = normrnd(mu1, sig1, 1, sum(complexAssignment == j));
					complexVolumes = max(4 / 3 * pi * complexRadii.^3, 0); % Disallow negative DNA/cell
					dnaPerCell(i, j) = sum(complexVolumes);
				end
			end
			
			% Scale to protein copies/AFU/MEFL/etc
			protPerCell = sf .* dnaPerCell;
			
			% Add autofluorescence
			protPerCell = protPerCell + normrnd(mu2, sig2, size(protPerCell));
		end
		
		
		function protPerCell = radiiRnd2(numCells, k, th, mu1, sig1, sf, mu2, sig2, numPoly)
			
			if ~exist('numPoly', 'var')
				numPoly = 1;
			end
			
			% Poisson for complexes per cell
			% See link below - sum of poissons = poisson w/ summed lambda
			% https://en.wikipedia.org/wiki/List_of_convolutions_of_probability_distributions
% 			complexesPerCell = nbinrnd(k/numPoly, 1/th, numCells, numPoly);
			complexesPerCell = round(gamrnd(k, th, numCells, 1));
			
			% Find amount of DNA (in terms of radii) per complex, then convert
			% to volume by relating radius to volume of a sphere
			dnaPerCell = zeros(numCells, numPoly);
			for i = 1:numCells
				complexAssignment = randi(numPoly, 1, complexesPerCell(i));
				for j = 1:numPoly
					complexRadii = normrnd(mu1, sig1, 1, sum(complexAssignment == j));
					complexVolumes = max(4 / 3 * pi * complexRadii.^3, 0); % Disallow negative DNA/cell
					dnaPerCell(i, j) = sum(complexVolumes);
				end
			end	
			
			% Scale to protein copies/AFU/MEFL/etc
			protPerCell = sf .* dnaPerCell;
			
			% Add autofluorescence
			protPerCell = protPerCell + normrnd(mu2, sig2, size(protPerCell));
		end
		
		
		function [PFit, minVal, pdfs] = fitPdfGammaNormal(xdata, P0, options)
			% Fits the PDF of data using convolved Gamma and Normal distributions
			%	
			%	[PFit, minVal, fitPdfs] = fitPdfGammaNormal(xdata, P0, options);
			%	
			%	PFit/P0 array order: ['k', 'th', 'mu', 'sigma', <'lambda'>]
			%						(lambda only used w/ 'poissGammaNorm' model)
			%
			%	'options' fields	(defaults | <options>)
			%		'model'			('gammaNorm' | 'gammaNorm', 'poissGammaNorm')
			%		'mode'			('log' | 'lin', 'log', 'both')
			%		'numRep'		(5)
			%		'numIter'		(1e3)
			%
			%	'pdfs' fields
			%		'xvalsLin'		Bin centers for linear PDF
			%		'dataLin'		PDF of data w/ linear binning
			%		'fitLin'		PDF of fit w/ linear binning
			%		'xvalsLog'		Bin centers for log PDF
			%		'dataLog'		PDF of data w/ log binning
			%		'fitLog'		PDF of fit w/ log binning
			%
			%	Note: density is logicle-transformed to ensure that the whole distribution 
			%	is weighted equally. *Give edgesLog in 'logicle' space*
			%
			%		Alt: apply a weight vector which is exponentially-decreasing
	
			% Check inputs
			zCheckInputs_fitPdfGammaNormal();
			
			% Find density of data in linear and logicle space
			data = struct();	% Global to avoid function overhead later
			makeSpaces();		% Adds centersLin/Log, edgesLin/Log, targetsLin
			data.densityLin = histcounts(xdata, data.edgesLin, ...
					'Normalization', 'Probability');
			data.densityLinShort = histcounts(xdata, data.edgesLinShort, ...
					'Normalization', 'Probability');
			[data.densityLog, data.validLog] = Fitting.adjustPdfBins( ...
					data.densityLin, data.centersLin, data.edgesLog);
			
			% Remove Inf and NaN values
			% --> We can't pre-emptively take these out of densityLin because
			%	  the spacing between points needs to be even for the fitting to
			%	  work due to the convolution
			data.validLin = (data.densityLin > 0);
			data.validLinShort = (data.densityLinShort > 0);
% 			data.validLog = (densityLog > 0); % Don't need - calculated where needed
			
% 			numel(data.densityLin)
% 			numel(data.centersLin)
% 			numel(data.validLin)

			% Find location of linear bins within logicle bins
			
			minVals = inf(1, options.numRep);
			P = repmat(P0, [1, options.numRep]);
			for ri = 1:options.numRep

				dP = ones(size(P0));

				for i = 1:options.numIter
					% dP is a change vector applied to the parameters
					%	It initizializes as 1, then changes according to whether 
					%	it makes the fit better or worse. Better = repeat change again
					%	Worse = try something else. 
					P0_i = P(:, ri) .* dP;
					fval = minfunc(P0_i);
					if (fval < minVals(ri))
						minVals(ri) = fval;
						P(:, ri) = P0_i;
					else
						dP(:) = 1;
						pi = randi(numel(P0));
						dP(pi) = 10.^(randn(1) / 5);
		% 				dP =  10.^(randn(size(dP)) / 5);
					end
				end
			end
		% 	minVals, P
			[minVal, minIdx] = min(minVals);
			PFit = P(:, minIdx);
			PFit_cell = num2cell(PFit);

			% Find parameters
		% 	options = optimset('Display', 'off');
		% 	[P, fval] = fminsearch(@minfunc, P0, options);

			% Generate final PDFs
			% --> Generate extra x-vals for the PDF which will be 
			%	  adjusted for log-space later
			switch options.model
				case 'gammaNorm'
					fitLin = Fitting.pdfGammaNormal( ...
							data.centersLin, PFit_cell{:});
					fitLinExtra = Fitting.pdfGammaNormal( ...
							data.centersLinExtra, PFit_cell{:});
				case 'poissGammaNorm'
					fitLin = Fitting.pdfPoissGammaNormal( ...
							data.centersLin, PFit_cell{:});
					fitLinExtra = Fitting.pdfPoissGammaNormal( ...
							data.centersLinExtra, PFit_cell{:});
			end	
			% Create log-space fit
% 			fitLog = Fitting.adjustPdfBins(fitLinExtra, ...
% 					data.centersLinExtra, data.edgesLog);
			fitLog = Fitting.adjustPdfBins(fitLin, ...
					data.centersLin, data.edgesLog);
			
			% Save linear outputs
			pdfs.xvalsLin = data.centersLin(data.validLin);
			pdfs.dataLin = data.densityLin(data.validLin);
			pdfs.fitLin = fitLin(data.validLin);
			
			% Save log outputs
			pdfs.xvalsLog = data.edgesLog(data.validLog);
			pdfs.dataLog = data.densityLog(data.validLog);
			pdfs.fitLog = fitLog(data.validLog);				
			
			fprintf(1, 'Fval = %.2f\n', minVal);


			% --- Helper Functions --- %


			function fval = minfunc(p)
				p_cell = num2cell(p);
				switch options.model
					case 'gammaNorm'
						fitFunc = @Fitting.pdfGammaNormal;
					case 'poissGammaNorm'
						fitFunc = @Fitting.pdfPoissGammaNormal;
				end	
		% 		xdiff = numel(xvals_lin) - numel(testPDF); % Needed when min(xvals) > 0 
				
				fval = 0;
				
				if ismember(options.mode, {'log', 'logicle', 'both'})
					testPdfLin = fitFunc(data.centersLin, p_cell{:});
					testPdfLog = Fitting.adjustPdfBins(testPdfLin, ...
							data.centersLin, data.edgesLog);	
					% --> Take validLog since we need it to adjust the data PDF
					fv = sum(abs(log10(testPdfLog(data.validLog)) - ...
							log10(data.densityLog(data.validLog))).^1);
					fval = fval + fv / sum(data.validLog); % Adjust for # points
				end
				if ismember(options.mode, {'lin', 'linear', 'both'})
					% The linear fit should happen on 'short' data (smaller vals), 
					% where the fit parameters matter more
					testPdfLinShort = fitFunc(data.centersLinShort, p_cell{:});
					fv = sum(abs(log10(testPdfLinShort(data.validLinShort)) - ...
							log10(data.densityLinShort(data.validLinShort))).^1);
					fval = fval + fv / sum(data.validLinShort); % Adjust for # points
					% Valid points in the log space are automatically determined, 
					% but the valid linear points need to be pulled in
				end
			end
			
			
			function makeSpaces()
				
				% Find range and step from data (all log-space)
				minSign = sign(min(xdata));
				expMin = ceil(log10(minSign * min(xdata)) * 10) / 10; % Force positive
				expMinPos = floor(log10(min(xdata(xdata > 0))) * 10) / 10;
				expMax = ceil(log10(max(xdata)) * 10) / 10;
				dx1 = 10^(expMin);
				dx2 = 10^(expMin - 1);
				
				% Find edges and centers of bins
				[data.centersLinShort, data.edgesLinShort] = Fitting.createSpace( ...
						minSign * 10.^expMin, 10.^(expMax - 1), dx2);
				
				[data.centersLin, data.edgesLin] = Fitting.createSpace( ...
						minSign * 10.^expMin, 10.^expMax, dx1);
				
				[data.centersLinExtra] = createSpace( ...
						minSign * 10.^expMin, 10.^expMax, dx2);
				
				[cL, eL] = Fitting.createSpace(expMinPos, expMax, 0.1); % 1/10 log steps
				data.centersLog = 10.^cL;
				data.edgesLog = 10.^eL;
			end
			

			function zCheckInputs_fitPdfGammaNormal()
				
				% X values
				validateattributes(xdata, {'numeric'}, {'vector'}, mfilename, 'xdata', 1);
				% Ensure xdata is non-inf and non-nan
				valid = (~isinf(xdata) & ~isnan(xdata));
				xdata = xdata(valid);
				
				% Model selection
				if isfield(options, 'model')
					validatestring(options.model, {'gammaNorm', 'poissGammaNorm'}, mfilename, 'model');
				else
					options.model = 'gammaNorm';
				end
				
				% Parameters
				switch options.model
					case 'gammaNorm'
						numParams = 4;
					case 'poissGammaNorm'
						numParams = 5;
				end
				validateattributes(P0, {'numeric'}, {'vector', 'numel', numParams}, mfilename, 'P0', 2);
				
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
					validatestring(options.mode, {'lin', 'linear', 'log', 'logicle', 'both'}, mfilename, 'mode');
				else
					options.mode = 'log';
				end				

			end
	
		end
		
	end
	
end