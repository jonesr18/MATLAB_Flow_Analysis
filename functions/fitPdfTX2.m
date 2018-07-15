function [lam, k, th, mu, sig, minVal, fitPDF] = fitPdfTX2(xvals_lin, xvals_log, density, gamOnly, numIter, doMEF, logicleParams)
	% Fits the transfection distribution to a sample x (including Poisson)
	%
	%	Note: density is compared on log-binned data to ensure that the whole
	%	distribution is weighted equally. 
	%
	%		Alt: apply a weight vector which is exponentially-decreasing
	
	% Check inputs
	gamOnly = exist('gamOnly', 'var') && gamOnly;
	if ~exist('numIter', 'var'), numIter = 3e3; end
	doMEF = (exist('doMEF', 'var') && doMEF);
	if ~exist('logicleParams', 'var'), logicleParams = struct(); end
	logicleParams.MEF = doMEF;
	
	[~, ~, bID] = histcounts(Transforms.lin2logicle(xvals_lin, doMEF, logicleParams), xvals_log);
	ubID = unique(bID(bID ~= 0));
	bIDs = cell(1, numel(ubID));
	for ubi = 1:numel(ubID)
		bIDs{ubi} = (bID == ubID(ubi));
	end
	
	% Find initial values for [lambda, k, theta, mu, sigma] by MC sim
	F0 = inf;
	P0 = [0.2; 0.05; 3e3; 0; 30]; 
	for j = 1:numIter / 10
		% RVs: lam: 0 - 1 | k: 0 - 1 | th: 1e2 - 1e4 | mu: -0.5 - 0.5 | sig: 1e1-1e2
		P0_j = [rand, rand, 10.^(2 * (rand + 1)), rand - 0.5, 10.^(rand + 1)];
		F0_j = minfunc(P0_j);
		if F0_j < F0
			P0 = P0_j;
			F0 = F0_j;
		end
	end	
	
	% Iterate to improve fit
	minVal = F0;
	dP = ones(size(P0));
	P = P0;
	for i = 1:numIter

		% dP is a change vector applied to the parameters
		%	It initizializes as 1, then changes according to whether 
		%	it makes the fit better or worse. Better = repeat change again
		%	Worse = try something else. 
		P0_i = P .* dP;
		fval = minfunc(P0_i);
		if (fval < minVal)
			minVal = fval; 
			P = P0_i;
		else
			if gamOnly, pi = randi(numel(P0) - 2); 
				else, pi = randi(numel(P0)); end
			dP(pi) =  10.^(randn(1) / 5);
		end
	end
	
	% Find parameters
% 	options = optimset('Display', 'off');
% 	[P, fval] = fminsearch(@minfunc, P0, options);
	
	% Unpack parameters
	Pc = num2cell(P);
	[lam, k, th, mu, sig] = Pc{:};
	
	% Generate final PDF
	fitPDF = pdfTX2(xvals_lin, lam, k, th, mu, sig);
	
	fprintf(1, 'Fval = %.2f\n', minVal);
	
	
	% --- Helper Functions --- %
	
	
	function fval = minfunc(p)
		testPDF = pdfTX(xvals_lin, p(1), p(2), p(3), p(4), p(5));
% 		xdiff = numel(xvals_lin) - numel(testPDF); % Needed when min(xvals) > 0 
		logTestPDF = makeLogPDF(testPDF);
		fval = sum(abs(log10(logTestPDF) - log10(density(ubID))).^2);
	end


	function logPDF = makeLogPDF(linPDF)
		logPDF = zeros(size(bIDs));
		for bi = 1:numel(bIDs)
			logPDF(bi) = sum(linPDF(bIDs{bi}));
		end
	end
end


