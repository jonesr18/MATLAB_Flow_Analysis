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
	%	Methods
	%
	%		getColormap(N)		Returns an Nx3 matrix of RGB colors
    %
    % Written By 
	% Ross Jones
    % jonesr18@mit.edu
    % Weiss Lab, MIT
    
    properties (Constant)
        
        COLORS = struct( ...
            'red',          struct('cm_top', [1.0, 0.6, 0.6], 'cm_mid', [1.0, 0.0, 0.0], 'cm_bot', [0.4, 0.0, 0.0]), ...
            'orange',       struct('cm_top', [1.0, 0.8, 0.6], 'cm_mid', [1.0, 0.5, 0.0], 'cm_bot', [0.4, 0.2, 0.0]), ...
            'yellow',       struct('cm_top', [1.0, 1.0, 0.4], 'cm_mid', [0.8, 0.8, 0.0], 'cm_bot', [0.2, 0.2, 0.0]), ...
            'chartreuse',   struct('cm_top', [0.8, 1.0, 0.6], 'cm_mid', [0.5, 1.0, 0.0], 'cm_bot', [0.2, 0.4, 0.0]), ...
            'green',        struct('cm_top', [0.6, 1.0, 0.6], 'cm_mid', [0.0, 1.0, 0.0], 'cm_bot', [0.0, 0.4, 0.0]), ...
            'teal',         struct('cm_top', [0.6, 1.0, 0.8], 'cm_mid', [0.0, 1.0, 0.5], 'cm_bot', [0.0, 0.4, 0.2]), ...
            'cyan',         struct('cm_top', [0.4, 1.0, 1.0], 'cm_mid', [0.0, 0.8, 0.8], 'cm_bot', [0.0, 0.2, 0.2]), ...
            'sky',          struct('cm_top', [0.6, 0.8, 1.0], 'cm_mid', [0.0, 0.5, 1.0], 'cm_bot', [0.0, 0.2, 0.4]), ...
            'blue',         struct('cm_top', [0.6, 0.6, 1.0], 'cm_mid', [0.0, 0.0, 1.0], 'cm_bot', [0.0, 0.0, 0.4]), ...
            'purple',       struct('cm_top', [0.8, 0.6, 1.0], 'cm_mid', [0.5, 0.0, 1.0], 'cm_bot', [0.2, 0.0, 0.4]), ...
            'magenta',      struct('cm_top', [1.0, 0.4, 1.0], 'cm_mid', [0.8, 0.0, 0.8], 'cm_bot', [0.2, 0.0, 0.2]), ...
            'hotpink',      struct('cm_top', [1.0, 0.6, 0.8], 'cm_mid', [1.0, 0.0, 0.5], 'cm_bot', [0.4, 0.0, 0.2]), ...
            'grey',         struct('cm_top', [0.8, 0.8, 0.8], 'cm_bot', [0.2, 0.2, 0.2]), ...
            'redblue',      struct('cm_top', [1.0, 0.0, 0.0], 'cm_mid', [0.85, 0.85, 0.85], 'cm_bot', [0.0, 0.5, 1.0]), ...
...%             'greenpurp',    struct('cm_top', [0.0, 1.0, 0.5], 'cm_mid', [0.85, 0.85, 0.85], 'cm_bot', [0.5, 0.0, 1.0]));
            'greenpurp',    struct('cm_top', [0.4, 1.0, 0.8], 'cm_mid', [0, 0, 0], 'cm_bot', [0.8, 0.4, 1.0]), ...
			'ocean',		struct('cm_top', [237 248 177] ./ 255, 'cm_mid', [65 182 196] ./ 255, 'cm_bot', [37 52 148] ./ 255));
     
    end
    
    
    properties (Access = public)
        
        matlab;
        cm_text;
        
    end
    
    
    methods (Access = public)
        
        function self = ColorMap(type)
            
            assert(ischar(type));
            
            if (exist(type) == 2) %#ok<EXIST>
                % colormap is MATLAB type
                self.matlab = true;
            else
                validatestring(type, fieldnames(self.COLORS), mfilename, 'type', 1);
                self.matlab = false;
            end
            
            self.cm_text = type;
        end
        
        
        function cm = getColormap(self, n, skew)
            % Returns a colormap with the specified length (n)
			%
			% 'skew' input allows for forcing the top/bottom values to
			% white/black: 'white', 'black', {'white', 'black'}
            
            % Check input
            validateattributes(n, {'numeric'}, {}, mfilename, 'n', 1);
            n = round(n);
            
            if (self.matlab)
                % If a MATLAB builtin, the colormap can be created directly using the name
                % of the colormap - here we use the eval function to get the colormap from
                % its name
                cm = eval(strcat(self.cm_text, '(', num2str(n), ')'));
            else
                % If not a MATLAB colormap, then we use one of our own
				top = self.COLORS.(self.cm_text).cm_top;
                bot = self.COLORS.(self.cm_text).cm_bot;
				if exist('skew', 'var')
					if ischar(skew), skew = {skew}; end
					if ismember(skew, 'white'), top = [1, 1, 1]; end
					if ismember(skew, 'black'), bot = [0, 0, 0]; end
				end
                
                % Setup colormap output
                cm = zeros(n, numel(top));
                
                % Check if the colormap has a middle value (eg white)
                if isfield(self.COLORS.(self.cm_text), 'cm_mid')
                    mid = self.COLORS.(self.cm_text).cm_mid;
                    
                    half_n = ceil(n / 2);
                    
                    for i = 1:numel(top)
                        cm(1:half_n, i) = linspace(bot(i), mid(i), half_n);
                        cm(half_n : n, i) = linspace(mid(i), top(i), n - half_n + 1);
                    end
                else
                    % Take care of normal case (no middle)
                    for i = 1:numel(top)
                        cm(:, i) = linspace(bot(i), top(i), n);
                    end
                end
            end
        end
        
    end
    
end