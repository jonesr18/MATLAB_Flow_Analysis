function out = gamma_ND(p, a)
	% Computes the multi-variate gamma function for p dimensions and input k
	
	out = pi^(p * (p - 1) / 4);
	for i = 1:p
		out = out * gamma(a + (1 - i) / 2);
	end
end