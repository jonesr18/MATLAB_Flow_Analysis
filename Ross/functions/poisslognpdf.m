function [plnpdf, ppdf, x] = poisslognpdf(x, lambda, mu, sigma, param)
	% Note: Inserts 0 into x if not present

	% Average x step value
	dx = mean(diff(x));
	
	% Insert 0 into x if necessary
	if ~any(x == 0)
		x1 = find(x > 0, 1);	% Only grab the first one
		x = [x(1:(x1 - 1)), 0, x(x1:end)];
	end
	x0 = (x == 0);
	
	% Get poiss pdf
	minProb = 1e-10;
	px = 0 : 1 : max(100, 20 * lambda);
	ppdf = poisspdf(px, lambda);
	pxint = find(ppdf > minProb);
	
	% Remove 0 
	pxint = pxint(2:numel(pxint));
	
	% Start combining logn * poiss w/ first integer
	%	Scale by probability of each integer to correctly weight
% 	plnpdf = ppdf(xint(1)) .* log10npdf(x, log10(x(xint(1))) + mu, sigma, param);
	
	% Convolve each additional logn * poiss to add subsequent complexes
	plnpdf = zeros(size(x));
	for xi = pxint
		
		W = ppdf(xi);			% Weight
		K = log10(px(xi));		% Scale factor
		plnpdf_i = W .* log10npdf(x, K + mu, sigma, param);
		
		plnpdf = plnpdf + plnpdf_i;
% 		plnpdf = fconv(plnpdf, plnpdf_i);
% 		plnpdf = plnpdf(1:numel(x)) .* dx; % Only take original x-vals
	end
	
	% Add 0 probability
	p0 = ppdf(1);
	pNon0 = 1 - p0;
	plnpdf = plnpdf ./ sum(abs(plnpdf)) ./ dx * pNon0;
	plnpdf(x0) = p0;
end


