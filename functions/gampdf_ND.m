function gpdfn = gampdf_ND(X, a, b, S)
	
% 	X = reshape(X, [], 1); 
	
	warning('Does not work currently');

	p = size(S, 1);
	assert(p == size(S, 2), 'C must be square!');
	
	gpdfn = det(S)^-a ./ (b^(p*a) * gamma_ND(p, a)) * ...
		det(X)^(a - (p + 1) / 2) * ...
		exp(trace(- 1/ b * S^-1 * X));
	
end