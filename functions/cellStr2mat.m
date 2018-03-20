function matrix = cellStr2mat(cellArray)
    % Converts a cell array of strings that are numbers to a MATLAB array of numbers.
    
    % Check input
    validateattributes(cellArray, {'cell'}, {}, mfilename, 'cellArray', 1);
    
    [H, W] = size(cellArray);
    matrix = zeros(H, W);
    for i = 1:H
        for j = 1:W
            matrix(i, j) = str2double(cellArray{i, j});
        end
    end
end