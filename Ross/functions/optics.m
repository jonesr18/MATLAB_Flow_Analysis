function [orderedPoints, reachDists] = optics(points, epsilon, minPts)
    % [ORDEREDPOINTS, REACHDISTS] = OPTICS(POINTS, EPSILON, MINPTS)
    % A simple OPTICS implementation of the original paper:
	% "OPTICS: Ordering Points To Identify the Clustering Structure (1999)
	% -- Mihael Ankerst, Markus M. Breunig, Hans-Peter Kriegel, Jörg Sander
	%
    % OPTICS handles differential cluster density better than DBSCAN
	%
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
    %
    %	Outputs: 
	%		
	%		orderedPoints:  A 1xP vector of indexes for points in the order they were processed
    %
    %       reachDists:     A 1xP vector of reachable distances for each point. 
    %                       Typically, this is plotted in order of orderedPoints
	%
	%
	% Written Ross Jones
    % jonesr18@mit.edu
	% github: jonesr18
    % Updated 2016-04-06 for clarity
	%
    % Credit to Tianxiao Jiang for original MATLAB DBSCAN implementation, 
    % which this is partially based off of.
    % 
    % Update log:
    
    %Initialize cluster membership as -1, which means UNCLASSIFIED
	checkInputs(points, epsilon, minPts)
	minPts = round(minPts); 		 		% Ensure integer
	P = size(points, 2);
	
	% Find distances
	distMat = dist(points);
	
	% Vector indicating if the point has been processed
    processed = false(1, P);
    
    % Vector of ordered points indicative of clusters
    orderedPoints = zeros(1, P);
    iii = 1;
    
    % Vector of reachability distances for each point
    reachDists = nan(1, P);
    
	% Randomly choose the visiting order and iterate over points
	visitOrder = randperm(P);
	for i = 1:P
		
		% For each point, check if it is unclassified before expanding
		pt = visitOrder(i);
		if ~processed(pt)
			% Point is not yet processed - mark that it is now being processed
            processed(pt) = true;
            
            % Add the point to the ordered list
            orderedPoints(iii) = pt;
            iii = iii + 1;
            
            % Find neighbors neaby pt and the core distance to minPts
            [neighbors, coreDist] = coreDistance(distMat(:, pt), epsilon, minPts);
            reachDists(pt) = coreDist;
            
            % If there are sufficient neighbors nearby, search them for further connections
			if ~isnan(coreDist)
				% Create new priority queue
                queue = PriorityQueue();
                
                % Update reachable distance in neighbors
                update(neighbors, coreDist, distMat(:, pt));
                
				% Iterate through elements in priority queue and check their neighbors for linkage
				while ~isempty(queue)
                    newPt = queue.dequeue();
                    [newNeighbors, newCoreDist] = coreDistance(distMat(:, newPt), epsilon, minPts);
                    processed(newPt) = true;
                    orderedPoints(iii) = newPt;
                    iii = iii + 1;
                    if ~isnan(newCoreDist)
                        update(newNeighbors, newCoreDist, distMat(:, newPt))
                    end
                end
                clear queue
			end
		end
    end
    
    
    % --- Helper Functions --- %
    
    
    function update(neighbors, coreDist, neighborhood)
        % Iterate over the given set of neighbors to update their reachable distance 
        % such that it is minimized.
        for j = reshape(neighbors, 1, []);
            if ~processed(j)
                newReachDist = max(coreDist, neighborhood(j));
                if isnan(reachDists(j))
                    % Reach distance for this point has not been defined yet
                    reachDists(j) = newReachDist;
                    queue.queue(newReachDist, j);
                elseif (newReachDist < reachDists(j))
                    reachDists(j) = newReachDist;
                    queue.dequeue(j);
                    queue.queue(newReachDist, j);
                end
            end
        end
    end


    function [neighbors, coreDist] = coreDistance(neighborhood, epsilon, minPts)
        % Finds all points within epsilon distance from a query point (its neighbors),
        % and returns their numerical indexes. 
        %
        % Also finds the core distance to the query point in the given neighborhood
        %
        % Neighborhood should be a column of distances of points to a query point, 
        % including the point itself.

        % Find neighbors
        neighbors = find(neighborhood <= epsilon);
        if (numel(neighbors) < minPts)
            coreDist = NaN;
        else
            coreDist = neighborhood(neighbors(minPts));
        end
    end


    function checkInputs(points, epsilon, minPts)
        % Checks inputs for validity

        validateattributes(points, {'numeric'}, {}, mfilename, 'points', 1);
        validateattributes(epsilon, {'numeric'}, {'scalar'}, mfilename, 'epsilon', 2);
        validateattributes(minPts, {'numeric'}, {'scalar'}, mfilename, 'minPts', 3);
    end
end


