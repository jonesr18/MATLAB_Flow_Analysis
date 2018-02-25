classdef CellList < handle
    % A class that implements some methods of the Java List interface using a
    % MATLAB cell array backbone
	%
	%   Implemented methods:
	%
	%		add			Add one item
	%		addAll		Add a list of items
	%		get			Get the value of an item at a given index
	%		remove		Remove an item from a given index
	%		set			Set the value of an item at a given index
	%
	%	MATLAB-specific methods:
	%		
	%		numel		Returns the number of items in the list
	%		length		Returns the length of the list
	%		toCell		Returns the items in the list as a cell array
	%		toVec		Returns the items in the list as a vector
	%		print		Prints the contents of the list
	%		randomize	Randomizes the location of the items in the list
	%
	% Written By
	% Ross Jones
	% jonesr18@mit.edu
	% Weiss Lab, MIT
	
    properties (SetAccess = private)
        data;
        size;
    end
    
    methods (Access = public)
        
        function self = CellList(capacity)
            
            % Check inputs
            if (exist('capacity', 'var'))
                validateattributes(capacity, {'numeric'}, {'positive'}, mfilename, 'capacity', 1);
                capacity = round(capacity);
            else
                capacity = 1000;
            end
            
            % Initialize cell structure which holds data
            self.data = cell(1, capacity);
            self.size = 0;
		end
        
		
        function numElements = numel(self)
            % Returns the current number of elements in the array - same as length!
            numElements = self.size;
		end
        
		
        function len = length(self)
            % Returns the current length of the array - same as numel!
            len = self.size;
        end
        
        
        function cells = toCell(self)
            % Returns the cell list as it is
            cells = self.data(1:self.size);
		end
		
		
		function vector = toVec(self)
			% Returns the cell list as a vector
			%
			%	NOTE: Only works if all data is numeric
			
			vector = cell2mat(self.toCell());
		end
        
        
        function self = add(self, item, index)
            % Adds the given item to the list.
            %   An index to insert at can be designated in the optional input.
            
            % Check inputs
            if (exist('index', 'var'))
                validateattributes(index, {'numeric'}, {'positive'}, mfilename, 'index', 2)
                index = round(index);
            else
                index = self.size + 1;
            end
            
            % Increase Capacity if necessary
            self.size = self.size + 1;
            if (self.size > numel(self.data))
                self.data = [self.data, cell(1, numel(self.data))];
            end
            self.checkIndex(index);
            
            % Handle special cases depending on where added (front, middle, back)
            if (index == 1)
                self.data = [{item}, self.data];
            elseif (index < self.size)
                self.data = [self.data(1:index - 1); {item}; self.data(index:end)];
            else % index == self.size
                self.data{index} = item;
            end
        end
        
        
        function self = addAll(self, itemArray, index)
            % Adds the given array of items to the list.
            %   An index to insert at can be designated in the optional input.
            %   Note that a single string will be treated as a character array 
            
            % Check inputs
            if (exist('index', 'var'))
                validateattributes(index, {'numeric'}, {'positive'}, mfilename, 'index', 2)
                index = round(index);
            else
                index = self.size + 1;
            end
            
            % Add all items
            for i = 1:length(itemArray)
                if iscell(itemArray)
                    self.add(itemArray{i}, index);
                else
                    self.add(itemArray(i), index);
                end
                index = index + 1;
            end
        end
        
        
        function item = get(self, index)
            % Returns the item at the end of the list.
            %   Optional input index returns the item at the given index.
            
            % Check input
            if (exist('index', 'var'))
                validateattributes(index, {'numeric'}, {'positive'}, mfilename, 'index', 1)
                index = round(index);
            else
                index = self.size;
            end
            self.checkIndex(index);
            
            item = self.data{index};
        end
        
        
        function item = remove(self, index)
            % Removes and returns the item at the end of the list.
            %   Optional input index returns the item at the given index.
            
            % Check input
            if (exist('index', 'var'))
                validateattributes(index, {'numeric'}, {'positive'}, mfilename, 'index', 1)
                index = round(index);
            else
                index = self.size;
            end
            self.checkIndex(index);
            
            % Fetch item
            item = self.data{index};
            
            % Physically remove item from list
            if (index == 1)
                self.data = self.data(2:end);
            elseif (index < self.size)
                self.data = [self.data(1:index - 1); self.data(index + 1:end)];
            else % index == self.size
                self.data{self.size} = {};
            end
            
            % Reduce size of cell array
            self.size = self.size - 1;
            if (self.size < numel(self.data) / 10)
                % If the number of elements becomes less than 1/10th the current capacity, 
                % reduce the size of the vector by 1/3 (making the current size ~1/3 full
                self.data = self.data(1:round(numel(self.data) / 3));
            end
        end
        
        
        function self = set(self, index, item)
            % Sets the given index to the given item.
            
            % Check inputs
            validateattributes(index, {'numeric'}, {'positive'}, mfilename, 'index', 1);
            self.checkIndex(index);
            
            self.data{index} = item;
        end
        
        
        function print(self)
            % Prints contents
            
            for i = 1:self.size
                disp(self.data{i})
            end
        end
        
        
        function self = randomize(self, seed)
            % Randomizes the current data in the list
            %   Optional input seed is for the random number generator
            
            % Check inputs
            if (exist('seed', 'var'))
                rng(seed)
            end
            
            shuffledIndexes = randperm(self.size);
            self.data = self.data(shuffledIndexes);
        end
        
    end
    
    
    methods (Access = private)
        
        function checkIndex(self, index)
            % Checks an index to ensure that it is valid.
            if (index > self.size)
                error(['Index too large, must be smaller than array size.\n' ...
                       'Size: %d\nIndex: %d\n'], self.size, index);
            end
        end
        
    end
end