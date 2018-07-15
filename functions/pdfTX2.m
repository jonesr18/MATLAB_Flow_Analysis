function [gapdf, ppdf, gpdf, apdf] = pdfTX2(x, lambda, k, theta, mu, sigma)
	
	assert(numel(unique(diff(x))) == 1, 'x must be evenly spaced for convolution')
	
	xpoiss = 0:100;
	ppdf = poisspdf(xpoiss, lambda); 
	gpdf = zeros(size(x)); 
	for pi = 1:numel(ppdf)
		gpdf = gpdf + ppdf(pi) * gampdf(x, k * xpoiss(pi), theta);
	end
	apdf = normpdf(x, mu, sigma);
	
	% Find x step value and generate x vals for conv vector
	dx = mean(diff(x));
	[~, ~, ~, x2_idx] = createSpace(min(x), max(x), dx);
	
	% Gamma pdf = inf at x = 0, so take that out and replace it w/ 
	% remaining probability from the summation of the rest
	gpdf = adjPdf(gpdf, dx);
	
	% Convolve and extract the conserved x values
	% --> Convolution basically extends the x vector 2*x in both directions from
	% 0, but since there is no x vector given, we have to adjust x manually.
	gapdf = fconv(gpdf, apdf);
	gapdf = gapdf(x2_idx) .* dx;
	
	% Normalize
	gapdf = gapdf ./ sum(gapdf);
end