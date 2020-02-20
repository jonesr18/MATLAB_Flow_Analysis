function output = collapse(matrix, dims)
	% Collapses the given matrix along the dimensions given by dims
	% 
	%	output = collapse(matrix, dims);
	%
	%	The new, collapsed dimension is placed in dimension min(dims)
	%
	%	The order of 'dims' will affect the order in the new dimension. 
	%
	% Written By
	% Ross Jones
	% jonesr18@mit.edu
	% Weiss Lab, MIT

	nonDims = setdiff(1:ndims(matrix), dims);
	matrixPermd = permute(matrix, [nonDims, dims]);
	
	sizes = cell(1, numel(nonDims));
	for si = 1:numel(sizes)
		sizes{si} = size(matrix, nonDims(si));
	end
	matrixReshd = reshape(matrixPermd, sizes{:}, []);
	
	dimsMR = ndims(matrixReshd);
	if (min(dims) == 1)
		% Handle case where one of the collapsed dimensions was the first
		order = [dimsMR, 1:dimsMR - 1];
	else
		% Handle all other cases
		order = [1:min(dims) - 1, dimsMR, min(dims):dimsMR - 1];
	end
	output = permute(matrixReshd, order);
end