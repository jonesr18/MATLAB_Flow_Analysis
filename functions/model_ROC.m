function truePosRate = model_ROC(p, falsePosRate)
	% This function is the model we use for fitting false pos rate vs true pos rate 
	% (1-specificity vs sensitivity) for ROC curves
	
	% Extract mean/std guess for positive population
	% Mean/std are for the evaluation metric - not necessarily the raw values
	mu2 = p(1);
	sigma2 = p(2);
	
	% Mean/std of negative population are assumed to be 0 +/- 1
	mu1 = 0;
	sigma1 = 1;
	
	probDist1 = makedist('Normal',mu1,sigma1);
	split = icdf(probDist1, 1 - falsePosRate);
	
	probDist2 = makedist('Normal',mu2,sigma2);
	truePosRate = 1 - cdf(probDist2, split);

end
