function [gapdf, gpdf, apdf] = pdfTX(x, k, theta, mu, sigma, C)
	
	assert(numel(unique(diff(x))) == 1, 'x must be evenly spaced for convolution')
	if exist('C', 'var')
		assert(size(C, 1) == size(C, 2), 'C must be square!'); 
	else
		C = 1;
	end
	
	gpdf = gampdf(x, k, theta);
% 	gpdf = gampdf_ND(x, k, theta, C);
	apdf = normpdf(x, mu, sigma);	% No covariance in autofluorescence
	
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
end