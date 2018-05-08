function [centers, edges, centers2, c2_idx] = createSpace(xmin, xmax, dx)
	% Generates centers, edges (for binning), the edges that result from
	% convolution of two centers-based pdf vectors, and the indexes
	% corresponding with centers in centers2. 

	edges = (xmin - dx / 2) : dx : (xmax + dx / 2);
	centers = xmin : dx : xmax;
	centers2 = ((2 * xmin) : dx : (2 * xmax)) - sign(xmin) * mod(xmin, dx);
	c2_idx = ismember(centers2, centers);
end