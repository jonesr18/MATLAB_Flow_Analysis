classdef ColorMap < handle
    % A class for building custom or MATLAB-designed colormaps
    %
	%	Can be given any MATLAB colormap or one of several novel colormaps:
    %   
    %		red
    %		orange
    %		yellow
    %		chartreuse
    %		green
    %		teal
    %		cyan
    %		sky
    %		blue
    %		purple
    %		magenta
    %		hotpink
    %		grey
	%		redblue
	%		greenpurp
	%
	%	matplotlib colormaps (From Stephen Cobeldick on MFEX)
	%
	%		viridis
	%		cividis
	%		inferno
	%		magma
	%		plasma
	%		twilight
	%		tab10
	%		tab20
	%		tab20b
	%		tab20c
	%
	%	Methods
	%
	%		getColormap(N)		Returns an Nx3 matrix of RGB colors
    %
    % Written By 
	% Ross Jones
    % jonesr18@mit.edu
    % Weiss Lab, MIT
    
    properties (Constant)
        
        HSV = struct( ...
            'red',		0   ./ 360, ...
            'orange',	30  ./ 360, ...
            'yellow',	60  ./ 360, ...
            'leaf',		90  ./ 360, ...
            'green',	120 ./ 360, ...
            'teal',		150 ./ 360, ...
            'cyan',		180 ./ 360, ...
            'sky',		210 ./ 360, ...
            'blue',		240 ./ 360, ...
            'purple',	270 ./ 360, ...
            'magenta',	300 ./ 360, ...
            'hotpink',	330 ./ 360);
		
		CUSTOM = struct( ...
            'redblue',      struct('top', [1.0, 0.0, 0.0], 'mid', [0.85, 0.85, 0.85], 'bot', [0.0, 0.0, 1.0]), ...
...%             'greenpurp',    struct('top', [0.0, 1.0, 0.5], 'mid', [0.85, 0.85, 0.85], 'bot', [0.5, 0.0, 1.0]));
            'greenpurp',    struct('top', [0.4, 1.0, 0.8], 'mid', [0, 0, 0], 'bot', [0.8, 0.4, 1.0]), ...				% Middle = black
			'greenpurp2',	struct('top', [0.0, 0.27, 0.11], 'mid', [0.85, 0.85, 0.85], 'bot', [0.17, 0.0, 0.21]), ...	% Middle = grey, from https://observablehq.com/@d3/color-schemes: PRGn
			'greenpurp3',	struct('top', [0.1, 0.6, 0.3], 'mid', [0.85, 0.85, 0.85], 'bot', [0.3, 0.1, 0.6]), ...		% Middle = grey, more vibrant colors
			'ocean',		struct('top', [237 248 177] ./ 255, 'mid', [65 182 196] ./ 255, 'bot', [37 52 148] ./ 255), ...
			'redK',			struct('top', [1.0, 0.0, 0.0], 'bot', [0.2, 0.0, 0.0]), ...
			'redW',			struct('top', [1.0, 0.0, 0.0], 'bot', [1.0, 0.8, 0.8]), ...
			'grnK',			struct('top', [0.5, 1.0, 0.0], 'bot', [0.1, 0.2, 0.0]), ...
			'grnW',			struct('top', [0.5, 1.0, 0.0], 'bot', [0.9, 1.0, 0.8]), ...
			'bluK',			struct('top', [0.0, 0.5, 1.0], 'bot', [0.0, 0.1, 0.2]), ...
			'bluW',			struct('top', [0.0, 0.5, 1.0], 'bot', [0.8, 0.9, 1.0]));
     
    end
    
    
    properties (Access = public)
        
        category = '';
        cmapName = '';
        
    end
    
    
    methods (Access = public)
        
        function self = ColorMap(cmapName)
            
            assert(ischar(cmapName));
            
            if (exist(cmapName) == 2) %#ok<EXIST>
                % colormap is a pre-built function (e.g. MATLAB built-in)
                self.category = 'function';
			elseif ismember(cmapName, fieldnames(self.HSV))
				self.category = 'hsv';
			elseif ismember(cmapName, fieldnames(self.CUSTOM))
				self.category = 'custom';
			else
				error('Colormap not recognized: %s', cmapName)
            end
            
            self.cmapName = cmapName;
        end
        
        
        function cm = getColormap(self, numColors, skew)
            % Returns a colormap with the specified length (n)
			%
			% 'skew' input allows for forcing the top/bottom values to
			% white/black: 'white', 'black', {'white', 'black'}
            
            % Check input
            validateattributes(numColors, {'numeric'}, {}, mfilename, 'n', 1);
            numColors = round(numColors);
            
            switch self.category
				case 'function'
					% If a MATLAB builtin, the colormap can be created directly using the name
					% of the colormap - here we use the eval function to get the colormap from
					% its name
					cm = eval(strcat(self.cmapName, '(', num2str(numColors), ')'));
				
				case 'hsv'
					% If HSV-based, we use the hue and just tune the sat/value
					hues = self.HSV.(self.cmapName) .* ones(numColors, 1);
					
					% A cleaner way is to use HSL, then convert to HSV, then to RGB
					lightnessMax = 0.8;
					lightnessMin = 0.2;
					
					if exist('skew', 'var')
						if ischar(skew), skew = {skew}; end
						if ismember('white', skew), lightnessMax = 1; end
						if ismember('black', skew), lightnessMin = 0; end
					end
					
					lightnesses = linspace(lightnessMin, lightnessMax, numColors)';
					
					% Convert from HSL to HSV
					%	L = 0.0 --> V = 0.0 | S = 1.0
					%	L = 0.5 --> V = 1.0 | S = 1.0
					%	L = 1.0 --> V = 1.0 | S = 0.0
					
					saturations = 1 - max((lightnesses - 0.5) * 2, 0);
					values = min(lightnesses * 2, 1);
					
					cm = hsv2rgb([hues, saturations, values]); 
				
				case 'custom'
					% If custom, we use a top/middle/bottom description scheme
					
					colorScheme = self.CUSTOM.(self.cmapName);
					top = colorScheme.top;
					bot = colorScheme.bot;

					% Setup colormap output
					cm = zeros(numColors, numel(top));

					% Check if the colormap has a middle value (eg white)
					if isfield(colorScheme, 'mid')
						mid = colorScheme.mid;
						
						if (numColors > 3)
							% Calculate a midpoint that maximally equalizes the gradient
							% of color change above and below the midpoint
							midDistBot = abs(mean(mid) - mean(bot));
							midDistTop = abs(mean(mid) - mean(top));
							optiMid = midDistBot / (midDistBot + midDistTop);

							midPoint = ceil(optiMid * numColors);
						elseif (numColors == 3)
							% If 3 colors, midpoint is middle
							midPoint = 2;
						else 
							% If 2 colors, midpoint is top, and if 1 color, just use midpoint. 
							midPoint = numColors;
						end

						for i = 1:numel(top)
							cm(1:midPoint, i) = linspace(bot(i), mid(i), midPoint);
							cm(midPoint : numColors, i) = linspace(mid(i), top(i), numColors - midPoint + 1);
						end
					else
						% Take care of normal case (no middle)
						for i = 1:numel(top)
							cm(:, i) = linspace(bot(i), top(i), numColors);
						end
					end
            end
        end
        
    end
    
end