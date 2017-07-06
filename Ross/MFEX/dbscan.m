function clust  =  dbscan(points, epsilon, minPts)
    % DBSCAN(POINTS, EPSILON, MINPTS)
	% A simple DBSCAN implementation of the original paper:
	% "A Density-Based Algorithm for Discovering clusters in Large Spatial
	% Databases with Noise" -- Martin Ester et.al.
	% Since no spatial access method is implemented,  the run time complexity
	% will be N^2 rather than N*logN
	%
	%**************************************************************************
	%	Inputs:
	%
	%		points:	 A DxP vector of points to cluster
    %                (D dimensions and P points)
	%
	%		epsilon: A scalar value for epsilonilon-neighborhood threshold.
	%
	%		minPts:  A scalar value for minimum points in 
	%				 epsilon-neighborhood that holds
	%				 the core-point condition.
	%
	%**************************************************************************
	%	Outputs: 
	%		
	%		clust:  An 1xP vector with cluster membership for each point. 
	%				0 is reserved for NOISE.
	%
	%**************************************************************************
	% Written by Tianxiao Jiang,  jtxinnocence@gmail.com
	% github: captainjtx
	% Nov-4-2015
	%
	% Updated 2016-04-06 for clarity
	% Ross Jones, jonesr18@mit.edu
	% github: jonesr18
	%
	%**************************************************************************

	%Initialize cluster membership as -1,  which means UNCLASSIFIED
	checkInputs(points, epsilon, minPts)
	minPts = round(minPts); 		 		% Ensure integer
	P = size(points, 2);
	
	% Find distances
	distMat = dist(points);
	
	% Set up cluster membership vector
	clust = repmat(-1, [1, P]);
	clusterId = 0;

	% Randomly choose the visiting order and iterate over points
	visitOrder = randperm(P);
	for i = 1:P
		
		% For each point, check if it is unclassified before expanding
		pt = visitOrder(i);
		if ~isclustered(clust(pt))
			% Point is not yet part of a cluster
			 
			neighbors = regionQuery(distMat(:, pt), epsilon);
			if 	(numel(neighbors) < minPts)
				% If the new cluster is noise, the cluster ID is set to zero so
				% that it is not part of any real cluster
				clust(pt) = 0;
			else
				% Increase cluster ID
				clusterId = clusterId + 1;
				
				% Add neighbors to cluster
				clust(neighbors) = clusterId;
				
				% Remove current point before expanding cluster
				neighbors = setxor(neighbors, pt);
				
				%Iteratively expand the cluster through density-reachability
				clust = expandCluster(clust, neighbors, clusterId, distMat, epsilon, minPts);
			end
		end
	end
end


function clust = expandCluster(clust, seeds, clusterId, distMat, epsilon, minPts)
	% Expand the cluster to find nearby neighbors to the point pt
	
	while ~isempty(seeds)
		
		% The currentP is removed at the end of the loop, so we can always
		% grab the next one off the top
		currentP = seeds(1);
		
		% Find neighbors to current point
		newSeeds = regionQuery(distMat(:, currentP), epsilon);
		
		% If there are enough neighbors, add to the cluster
		if numel(newSeeds) >= minPts
			for i = 1:numel(newSeeds)
				newP = newSeeds(i);
				if ~isclustered(clust(newP))
					if clust(newP) == -1 % not expanded yet
						seeds = [seeds(:) ; newP];
					end
					clust(newP) = clusterId;
				end
			end
		end
		seeds = setxor(seeds, currentP);
	end
end


function clustered = isclustered(queryPoint)
	% Returns true if the point
	clustered = (queryPoint > 0);
end


function neighbors = regionQuery(neighborhood, epsilon)
	% Finds all points within epsilon distance from pt, and returns their
	% numerical indexes 
	
	% Perform region query
	neighbors = find(neighborhood <= epsilon);	
end


function checkInputs(points, epsilon, minPts)
	% Checks inputs for validity
	
	validateattributes(points, {'numeric'}, {}, mfilename, 'points', 1);
	validateattributes(epsilon, {'numeric'}, {'scalar'}, mfilename, 'epsilon', 2);
	validateattributes(minPts, {'numeric'}, {'scalar'}, mfilename, 'minPts', 3);
end