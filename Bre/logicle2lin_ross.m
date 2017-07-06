function out=logicle2lin(X)
% LOGICLE2LIN(X) converts a vector of values from a logicle scale to a linear scale 
%   OUT = LOGICLE2LIN(X) returns a vector of values corresponding to the values 
%   in X transformed from a logicle scale to a linear one using
%   the conversion function described by Parks, et al. "A New Logicle
%   Display Method Avoids Deceptive Effects of Logarithmic Scaling for Low
%   Signals and Compensated Data" Equation (5), where
%     X - a vector of logicle values
%     out - a vector of linear converted values from X
%
%   Example:
%       out = logicle2lin(linspace(0,4));
%
%   Written by
%   Breanna Stillo
%   bstillo@mit.edu
%   Last Updated: 2014-10-14
%
%   Edit 2.9.15 by Ross Jones
%   Added comments and intuitive variable names
%   Major speed up by using algebraic conversion instead of symbolic toolbox

    % Set values
    DATA_MAX = 2^18;         % (T) The maximum measurable data value ( = 26214 )
    DISPLAY_BREADTH = 4.5;   % (M) Breadth of the display in decades
    DATA_MIN = -1500;        % (r) Negative range reference value

    % This gives a range for linearization around zero (W).
    linRange = (DISPLAY_BREADTH - log10(DATA_MAX / abs(DATA_MIN))) / 2;

    % p and W (linRange) are considered one parameter, p is introduced by the 
    % authors for compactness. Here we find p that solves their equivalency:
    %   w = 2p ln(p)/(p + 1)
    % Solved by WolframAlpha:
    P = linRange / (2 * lambertw(0.5 * exp(-linRange / 2) * linRange));

    % Find where the given data is out of the linear range.
    posX = (X >= linRange);

    % Compute and return the final linearized vector 
    out = toLin(X) .* posX - toLin(2 * linRange - X) .* (1 - posX);
    
    % This inner function does a final conversion
    function converted = toLin(convert)
        converted = DATA_MAX .* 10^-(DISPLAY_BREADTH - linRange) .* ...
            (10.^(convert - linRange) - P^2 .* 10.^(-(convert - linRange) ./ P) + P^2 - 1);
    end
end




