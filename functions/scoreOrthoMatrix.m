function [score, maxVals] = scoreOrthoMatrix(subMatrix)
	% From Bre:
	%  "This is the function that evaluates a certain square submatrix… 
	%	trying to figure out what the “score” is.  
	%	
	%	Looks like I took the value of a current “ON” sample and divided 
	%	by the values in it’s respective row and column and then multiplied 
	%	those together.  
	%
	%	Then multiplied those values for each of the “ON” samples together.  
	%
	%	I had tried a bunch of things and this seemed to work the best 
	%	for colors.

	maxVals = false(size(subMatrix));	% Store 1 where max values found
	M = subMatrix;						% Updates over time to set non-max 
											% values in the same row/col to NaN
	
	% Find r & c of max val in matrix, set the rest of row and col to
	% NaN and set that index to 1 in maxVals
	%	This seems to try finding the strongest val, then the next, and so
	%	on while eliminating concurrent rows/vals
	while (sum(maxVals(:)) < size(subMatrix, 1))
		[rowMax, colMax] = ind2sub(size(M), find(M == max(M(:))));
		maxVals(rowMax, colMax) = 1;
		M(rowMax, :) = NaN;
		M(:, colMax) = NaN;
	end 

	% Normalize each max value identified above by the other values in its 
	% row and column to identify those that are the most "unique". 
	maxValsIdxs = find(maxVals);
	maxValsScores = zeros(1, length(maxVals));
	for i = 1:length(maxValsIdxs)
		
		% Extract a max value
		[rowIdx, colIdx] = ind2sub(size(maxVals), maxValsIdxs(i));
		
		% Normalize to row and column it is in
		rowDiffs = subMatrix(rowIdx, colIdx) ./ subMatrix(rowIdx, :);
		colDiffs = subMatrix(rowIdx, colIdx) ./ subMatrix(:, colIdx);
		
		% Fail the submatrix if any max value is not the max in its row/col
		if any((rowDiffs < 1) | (colDiffs < 1))
			score = NaN;
			return
		else
			% Otherwise, score the submatrix max val-by-max val. 
			maxValsScores(i) = prod(rowDiffs(:)) * prod(colDiffs(:));
		end
	end
	
	% Final score
	score = prod(maxValsScores);
end




