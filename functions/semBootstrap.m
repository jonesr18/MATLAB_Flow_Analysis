function semb = semBootstrap(data, iterations)
	% Computes the standard error of the mean on a single sample using a
	% bootstrapping method with the given number of iterations (default = 1000);
	%
	%	Inputs
	%		data		<numeric> Numerical vector of data to calculate SEM on
	%		iterations	(optional) <numeric> The number of bootstrap iterations to run
	%
	%	Outputs
	%		sem			<numeric> The computed SEM. Returns NaN for empty inputs.
	
	% Return NaN for empty inputs
	if isempty(data), semb = NaN; return; end
	
	% Check inputs
	validateattributes(data, {'numeric'}, {'vector'}, mfilename, 'data', 1);
	if exist('iterations', 'var')
		validateattributes(iterations, {'numeric'}, {'scalar', 'positive'}, mfilename, 'iterations', 2);
		iterations = round(iterations);
	else
		iterations = 1000;
	end
	
	% Do bootstrapped mean calculation
	means = zeros(1, numel(iterations));
	for i = 1:iterations
		% Subsample
		idxs = randi(numel(data), 1, numel(data));
		
		% Compute mean
		means(i) = mean(data(idxs));
	end
	
	% Compute SEM 
	semb = std(means);
end