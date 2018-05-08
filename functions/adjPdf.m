function xpdf = adjPdf(xpdf, dx)
	% Replaces inf values in a pdf with the value that should be there, assuming
	% the entire pdf * dx sums to 1
	
	x0 = isinf(xpdf);
	p0 = (1 - sum(xpdf(~x0)) .* dx) ./ dx;	% Adjust so sum * dx = 1
	xpdf(x0) = p0 ./ sum(x0);
	
end