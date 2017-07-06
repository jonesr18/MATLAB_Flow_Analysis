%  Priority Queue class for Matlab
%   
%   Methods Implemented:
%      
%       peek        Look at the top-priority element w/o removing
%
%       dequeue     Remove the top-priority element
%
%       queue       Add a new element to the queue
%
%       remove      Remove an element (based on matching user data)
%
%   This class uses a binary heap implemented in a MATLAB array to prioritize data
%
%
% Written by Richard T. Guy
%
% Edited by Ross Jones to include removal of non-top priority nodes
% Also added more/better comments
classdef PriorityQueue < handle
   
    properties (SetAccess = private, GetAccess = public)
        size; 
        heapData = [];
        userData = {};  
    end
    
    methods
       
        function obj = PriorityQueue(capacity)
            % Initialize the priority queue with the given (optional) capacity. 
            % The capacity defaults to 1000 and is doubled every time it reaches capacity
            % and reduced 1/3 every time it drops to 10%
            
            if ~exist('capacity', 'var')
                capacity = 1000;
            end
            
            obj.size = 0;
            obj.userData = cell(capacity, 1);
            obj.heapData = zeros(capacity, 1);
        end
       
        
        function emptiness = isempty(obj)
            emptiness = (obj.size == 0);
        end
        
        
        function s = getSize(obj)
            s = obj.size; 
        end
                
        
        function [element_data, element_priority] = peek(obj)
            % Peek at the top-priority element w/o removing
            
            if (obj.size == 0)
                error('Binary heap error: accessed empty heap.');
            end  
            
            element_data = obj.userData{1};
            element_priority = obj.heapData(1);
        end
        
        
        function [element_data, element_priority] = dequeue(obj, element_data)
            % Dequeue the best element
            %
            %   User can also optionally supply an element to find and pop 
            %   This is slow and O(N) rather than O(1) as pop should be, since
            %   the heap is not optimized for random access.
            
            if (obj.size == 0)
                error('Binary heap error: accessed empty heap.');
            end    
            
            % Check if a particular piece of data was requested, otherwise pop the 
            % top-priority element
            if exist('element_data', 'var')
                index = obj.findMember(element_data);
            else
                index = 1;
            end
            
            % Access return data
            element_data = obj.userData{index};
            element_priority = obj.heapData(index);

            % Now percolate down to find a hole to move the terminal value
            hole = obj.percolateDown(index, obj.heapData(obj.size));
            obj.heapData(hole) = obj.heapData(obj.size);
            obj.userData(hole) = obj.userData(obj.size);
            
            % Reduce size of cell array (if necessary)
            obj.size = obj.size - 1;
            if (obj.size < numel(obj.heapData) / 10)
                % If the number of elements becomes less than 1/10th the current capacity, 
                % reduce the size of the vector by 1/3 (making the current size ~1/3 full
                obj.userData = obj.userData(1:round(numel(obj.userData) / 3));
                obj.heapData = obj.heapData(1:round(numel(obj.heapData) / 3));
            end
        end
        
        
        function queue(obj, element_priority, element_data)
            % Queue element onto heap
            %
            % element_priority is the "score" that this heap uses.
            % element_data should be a UNIQUE identifier
            %   - otherwise pop will not work for requesting based on element_data
            %   - can be a string OR a number, either will work!
            
            % Increase Capacity if necessary
            obj.size = obj.size + 1;
            if (obj.size > numel(obj.heapData))
                obj.userData = [obj.userData; cell(obj.size, 1)];
                obj.heapData = [obj.heapData; zeros(obj.size, 1)];
            end
            
            % Percolate up to find where to insert data
            hole = obj.percolateUp(obj.size, element_priority);
            obj.heapData(hole) = element_priority;
            obj.userData{hole} = element_data;
        end 
      
        
        function clear(obj)
            % Clears the object properties (reset)
            
            obj.size = 0;
            obj.userData = cell(1000, 1);
            obj.heapData = zeros(1000, 1);
        end
        
    end
  
    methods (Access=private)
        
        function index = findMember(obj, element_data)
            % Finds the given string or number in userData, and returns the index
            
            for i = 1:obj.size
                if ischar(element_data)
                    if strcmpi(element_data, obj.userData{i})
                        index = i;
                        return
                    end
                elseif isnumeric(element_data)
                    if (element_data == obj.userData{i})
                        index = i;
                        return
                    end
                end
            end
            error('Element not found')
        end
        
        
        function hole = percolateUp(obj, hole, item)
            % Percolate up through the array to find the hole where a new value is 
            % to be inserted
            
            % Swap values and recurse if "parent" has a higher priority value  
            target = floor(hole / 2);
%             fprintf(1, 'item: %d\nhole: %d\ntarget: %d\n', item, hole, target);
            if (hole > 1 && item < obj.heapData(target))
                obj.heapData(hole) = obj.heapData(target);
                obj.userData(hole) = obj.userData(target);
                hole = obj.percolateUp(target, item);
            end
        end
            
        
        
        function hole = percolateDown(obj, hole, item)
            % Percolate down through the array to find the hole where the terminal 
            % value (item) is to be moved
            
            if (2 * hole <= obj.size)
                % Find right and left "children"
                left = 2 * hole;
                right = left + 1;
                
                % Set target for next recursion
                target = right;
                if (right > obj.size || obj.heapData(left) < obj.heapData(right))
                    target = left;
                end
                
                % Swap values and recurse if target has a lower priority value
                if (obj.heapData(target) < item)
                    obj.heapData(hole) = obj.heapData(target);
                    obj.userData(hole) = obj.userData(target);
                    hole = obj.percolateDown(target, item);
                end
            end
        end
    end
end