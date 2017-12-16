classdef Gating < handle
	% Compilation of gating procedures from the Weiss Lab flow cytometry MATLAB repository
	%
	%	Methods are implemented as static functions, so they can be called directly
	%	or after creating a handle to this class.
	%
	%		Ex:
	%		[gData, gIdx] = Gating.gate(...);
	%
	%		Ex2:
	%		g = Transforms();
	%		[cells, vert] = g.gatePolygon(x, y, axScales, pos)
	
	properties (Constant)
		SCATTER_CHANNELS = {'SSC_A', 'SSC_H', 'SSC_W', 'FSC_A', 'FSC_H', 'FSC_W'};
	end
	
	methods (Static)
	
		function [gateData, gateIndices] = gate(data, indices, channelX, channelY, plotScale)
            % Runs a figure to let the user create a gate for the given channels by drawing 
            % a lasso around the data points desired. All plots are log-log scaled
            %
            %   Inputs:
            %       data        Data structure
            %       indices     indices of the data points to gate on from the input channels
            %       channelX    The name of the channel on the x-axis
            %       channelY    The name of the channel on the y-axis
            %       plotScale   The type of plot to use ('plot', 'semilogx', 'semilogy', 'loglog', 'biexp')
            %
            %   Outputs:
            %       gateData    A struct with the following fields:
            %           xlasso      The x-points defining the gate polygon
            %           ylasso      The y-points defining the gate polygon
            %                       NOTE: These will be in logicle format if plotScale is 'biexp'
            %           parent      The name of the parent gate ('none' if root). 
            %                       NOTE: this function sets this field to 'none'
            %           pctParent   The percentage of the parent
            %       gateIndices A list of indices of points passing the gate
            %
            %   Written by
            %   Ross Jones
            %   jonesr18@mit.edu
            %   Updated 2015-02-09

            % Check inputs
            plotScale = lower(plotScale);
            checkInputs_gate(data, indices, channelX, channelY, plotScale);

            % For standardization and ease of use, convert to logical indexing
            if (~islogical(indices))
                idx = false(1, data(1).nObs);
                idx(indices) = true;
                indices = idx;
            end

            % Define axes parameters
            AXES_MIN = 1;
            AXES_MAX = 2^18;
            biexp = false;

            % Set plot axes scales
            figure()
            ax = gca;
            hold on
            switch (plotScale)
                case 'plot'
                    ax.XScale = 'linear';
                    ax.YScale = 'linear';
                case 'semilogx'
                    ax.XScale = 'log';
                    ax.YScale = 'linear';
                case 'semilogy'
                    ax.XScale = 'linear';
                    ax.YScale = 'log';
                case 'loglog'
                    ax.XScale = 'log';
                    ax.YScale = 'log';
                case 'biexp'
                    biexp = true;
                    ax.XScale = 'linear';
                    ax.YScale = 'linear';
                    AXES_MIN = lin2logicle(-150);
            end
            ax.XLim = [AXES_MIN, AXES_MAX];
            ax.YLim = [AXES_MIN, AXES_MAX];

            % Plot the channels given and let the user draw a lasso
            % around the points that they would like to use, defining a gate
            xdata = data.(channelX).raw(indices);
            ydata = data.(channelY).raw(indices);
            if (biexp)
                Plotting.biexplot(xdata(1:10:end), ydata(1:10:end), 'density', 1)
            else
                Plotting.densityplot(xdata(1:10:end), ydata(1:10:end));
            end
            xlabel(channelX)
            ylabel(channelY)
            title('Click and drag to draw a lasso around the points to gate')
            [gateIndices, ~, ~, xlasso, ylasso] = selectdata('SelectionMode', 'lasso');

            % gateIndices is given as numerical indices of the sub-population of data, so
            % in order to make these indices mean anything they must be converted to 
            % indices of the entire population.

            % To start, we get the numerical indexes of the sub-population
            popIndices = find(indices);

            % Then we find the sub-population that passes
            passedPop = popIndices(gateIndices);

            % Then, for simplicity outside of this function, we convert back to logical indexes
            gateIndices = false(1, data(1).nObs);
            gateIndices(passedPop) = true;

            % Save the polygon defining the gate and population statistics
            gateData = struct( ...
                'xlasso', [], ...
                'ylasso', [], ...
                'indices', [], ...
                'parent', 'none', ...
                'pctParent', []);
            [gateData.xlasso, gateData.ylasso] = connectLasso(xlasso, ylasso);
            gateData.pctParent = 100 * sum(gateIndices) / sum(indices);

            % Update figure to include gate
            passX = data.(channelX).raw(gateIndices);
            passY = data.(channelY).raw(gateIndices);
            hold on     % have to turn this back on because of selectdata()
            if (biexp)
                plot(lin2logicle(passX), lin2logicle(passY), '.')
            else
                plot(passX, passY, '.')
            end
            plot(gateData.xlasso, gateData.ylasso)
            title('Gated Data')
            hold off


            % --- Helper Functions --- %

            function [xlasso, ylasso] = connectLasso(xlasso, ylasso)
                % Uses the first and lass coordinates of the lasso to connect it back together.
                % The connection is made with 50 points, just to give some resolution.

                NUM_POINTS = 50;

                % Collect the first and last points
                xfirst = xlasso(1);
                yfirst = ylasso(1);
                xlast = xlasso(end);
                ylast = ylasso(end);

                % x values are simply evenly spaced from the first to last
                xnew = linspace(xfirst, xlast, NUM_POINTS);

                % y values are calculated based on the slope from the first to last x values
                slope = (ylast - yfirst) / (xlast - xfirst);
                ynew = yfirst + slope * (xnew - xfirst);

                % Add values to xlasso and ylasso
                xlasso = [xlasso, xnew];
                ylasso = [ylasso, ynew];
            end


            function checkInputs_gate(data, indices, channelX, channelY, plotScale)
                % Checks the inputs to the main function

                % data
                validateattributes( ...
                    data, {'struct'}, {}, ...
                    mfilename, 'data', 1);

                % indices
                validateattributes( ...
                    indices, {'numeric', 'logical'}, {'integer'}, ...
                    mfilename, 'indices', 2);

                % channelX
                validatestring( ...
                    channelX, data.chanNames, ...
                    mfilename, 'channelX', 3);

                % channelY
                validatestring( ...
                    channelY, data.chanNames, ...
                    mfilename, 'channelY', 4);

                % plotScale
                validatestring( ...
                    plotScale, {'plot', 'semilogy', 'semilogx', 'loglog', 'biexp'}, ...
                    mfilename, 'plotScale', 6);
            end
			
        end
		
		
        function [gateP1, gateP2, gateP3] = standardGating(data, onlyP1, swap, maxPoints)
            % Creates 3 standard gates from a given samples
            %
            %   Gates
            %   
            %       P1:     FSC_A vs SSC_A
            %       P2:     FSC_H vs FSC_W
            %       P3:     SSC_H vs SSC_W
            %
            %   Inputs
            %
            %       data		A standard struct (must contain channels noted above)
			%		onlyP1		(optional) Logical flag to only gate on P1. The other 
			%					gates are returned as empty so that the output 
			%					interface is conserved. <Default := false>
			%		swap		(optional) Logical flag to swap the FSC and SSC
			%					plot axes (useful for Koch LSRII-HTS2)
			%		maxPoints	(optional) Number indicating the maximum number of 
			%					points to plot for drawing gates. <Default := 20000>
            %
            %   Outputs
            %
            %       gateP1      Logical index array for events in P1
            %       gateP2      Logical index array for events in P2
            %       gateP3      Logical index array for events in P3
            %
            %
            % Written by Ross Jones
            % Weiss Lab, MIT
            % Update Log: 
			%	2017-09-25
			%		Added data combining so that all samples can be passed together
			
			% Check inputs
			validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
			missingChannels = setdiff(Gating.SCATTER_CHANNELS, fieldnames(data));
			assert(isempty(missingChannels), 'Scatter channel missing: %s\n', missingChannels{:});
			onlyP1 = (exist('onlyP1', 'var') && all(logical(onlyP1))); % Default FALSE
			swap = (exist('swap', 'var') && all(logical(swap)));
			if ~exist('maxPoints', 'var'), maxPoints = 20000; end
			
			% For each scatter channel, extract # cells/# samples from each sample
			combinedData = struct();
			for ch = Gating.SCATTER_CHANNELS
				combData = [];
				for d = 1:numel(data)
					if ~isempty(data(d).(ch{:})) % Some FlowData tcData will be empty 
						combData = [combData; data(d).(ch{:}).raw(1:numel(data):end)]; %#ok<AGROW>
					end
				end
				combinedData.(ch{:}).raw = combData;
			end
			totalPoints = numel(combinedData.(Gating.SCATTER_CHANNELS{1}).raw);
			numPoints = min(maxPoints, totalPoints);
            
            % Full-size index for whole set of observations
            % In order to avoid problems with sub-indexing below, we will copy the entire FALSE index array
            % when assigning data to a new gate index array, then add back the TRUE entries by numerical index
            falseIdx = false(numPoints, 1); 
			pointsIdx = randperm(totalPoints, numPoints);

            % First do FSC-A vs SSC-A
            [subIdxP1, gateP1] = Gating.gatePolygon( ...
					combinedData.FSC_A.raw(pointsIdx), ...
					combinedData.SSC_A.raw(pointsIdx), 'semilogy'); 
			
			if onlyP1
				gateP2 = [];
				gateP3 = [];
			else
				if swap
					chans = {'FSC_H', 'FSC_W', 'SSC_H', 'SSC_W'};
				else
					chans = {'FSC_W', 'FSC_H', 'SSC_W', 'SSC_H'};
				end
				
				idxP1 = subIdxP1;               % No sub-indexing for P1
				numIdxP1 = find(idxP1);         % Find numerical indexes for P1

				% Second do FSC-W vs FSC-H
				[subIdxP2, gateP2] = Gating.gatePolygon( ... 
						combinedData.(chans{1}).raw(idxP1), ...
						combinedData.(chans{2}).raw(idxP1), 'linear');
				
				numIdxP2 = numIdxP1(subIdxP2);  % Get positions in P1 of objects passing P2 gate
				idxP2 = falseIdx;               % Copy this to get the full-sized indexing vector
				idxP2(numIdxP2) = true;         % Assign TRUE values according to what passes P1 AND P2

				% Third do SSC-W vs SSC-H
				[~, gateP3] = Gating.gatePolygon( ... 
					combinedData.(chans{3}).raw(idxP2), ...
					combinedData.(chans{4}).raw(idxP2), 'semilogy');    
			end
        end
        
        
        function data = applyStandardGates(data, gateP1, gateP2, gateP3, swap)
            % Applies the standard gates (see Gating.standardGating()) to the given sample
            % 
			%	Inputs
			%
			%		data		Data in a standard struct
			%		gateP1, gateP2, gateP3	...	
			%					Gate polygon outputs of Gating.standardGating()
			%		swap		(optional) Logical flag to swap the FSC and SSC
			%					plot axes (useful for Koch LSRII-HTS2)
			%
			%	Outputs
			%		data		Updated data struct w/ gate logicals added in 
			%					'gates' field.
            
			swap = (exist('swap', 'var') && all(logical(swap)));
			
            % Iterate over all samples in struct
            for d = 1:numel(data)
                
				if isempty(data(d).FSC_A) % Some FlowData tcData will be empty 
					continue
				end
				
                % Full-size false vector for setting up gates
                falseIdx = false(data(d).nObs, 1);

                % Apply gates
                [subIdxP1, ~] = Gating.gatePolygon( ...
                    data(d).FSC_A.raw, ...
                    data(d).SSC_A.raw, ...
                    'semilogy', gateP1);
                numIdxP1 = find(subIdxP1);
                idxP1 = falseIdx;
                idxP1(numIdxP1) = true;
				
				% Save Gate P1
				data(d).gates.P1 = idxP1;
				
				if swap
					chans = {'FSC_H', 'FSC_W', 'SSC_H', 'SSC_W'};
				else
					chans = {'FSC_W', 'FSC_H', 'SSC_W', 'SSC_H'};
				end
				
				% Gates P2/3 are not always applied
				if (exist('gateP2', 'var') && ~isempty(gateP2))
					[subIdxP2, ~] = Gating.gatePolygon( ...
						data(d).(chans{1}).raw(idxP1), ...
						data(d).(chans{2}).raw(idxP1), ...
						'linear', gateP2);
					numIdxP2 = numIdxP1(subIdxP2);
					idxP2 = falseIdx;
					idxP2(numIdxP2) = true;
					
					% Save Gate P2
					data(d).gates.P2 = idxP2;
				end
				
				if (exist('gateP3', 'var') && ~isempty(gateP3))
					[subIdxP3, ~] = Gating.gatePolygon( ...
						data(d).(chans{3}).raw(idxP2), ...
						data(d).(chans{4}).raw(idxP2), ...
						'semilogy', gateP3);
					numIdxP3 = numIdxP2(subIdxP3);
					idxP3 = falseIdx;
					idxP3(numIdxP3) = true;
					
					% Save Gate P3
					data(d).gates.P3 = idxP3;
				end
            end
        end
        
        
		function [cellsWithinGate, gateVertices] = gatePolygon(xdata, ydata, axis_scales, position)
			% gates a given population.  Uses the polygon vertices specified by
			% 'position' or if this is not given, the user is asked to specify the polygon
			% through UI
            
            % axes limits
            AXES_MIN = 1;
            AXES_MAX = 2^18;
                        
            % position input is optional, and thus the method will only pull up a GUI to create a 
            % new gate if it is not included
            
			if ~exist('position', 'var')
			    figure()
                ax = gca();
				if strcmp(axis_scales,'loglog')
					loglog(xdata,ydata,'.','MarkerSize',2)
					ax.XScale = 'log';
                    ax.YScale = 'log';
				elseif strcmp(axis_scales,'semilogy')
					%dscatter(xdata,ydata,'MARKER','.','MSIZE',2)
					semilogy(xdata,ydata,'.','MarkerSize',2)
					ax.YScale = 'log';
				elseif strcmp(axis_scales,'semilogx')
					semilogx(xdata,ydata,'.','MarkerSize',2)
					ax.XScale = 'log';
				else
					plot(xdata,ydata,'.','MarkerSize',2)
				end
				hold(ax, 'on')
                ax.XLim = [AXES_MIN, AXES_MAX];
                ax.YLim = [AXES_MIN, AXES_MAX];
                
				h = impoly(ax, 'Closed', true); % modified by JG to use impoly to determine if event is within gate - much easier this way.
				position = wait(h);
			end

			gateVertices = position;


			% determine if in or out of polygon
			cellsWithinGate = inpolygon(xdata, ydata, position(:, 1), position(:, 2));
            
		end
		
		
		function fcsdatGated = applyJCGate(sampleFile, gateFile, time_indicator) %modified by JG for time gating. (see applyJCGate for more comments)
			% applies a certain JC gate to a specified sample

			if nargin >= 3 && ~isempty(strfind(time_indicator,'yes')) %modified by JG only use time gating if user inputs 'yes' as 3rd input
				use_time = 1;
			else
				use_time = 0;
			end

			suffix = gateFile(end-2 : end);
			if strcmp(suffix, 'fcs')
				Gating.createJCGate(gateFile, use_time); %createJCGate checks to see if there is a mat file and doesn't run if there is one
				baseName = gateFile(1:end-4);
				gateFileName = [baseName '.mat'];% modified by JG. previously was: gateFileName=['JCGate_' baseName '.mat'];

				if exist(['JCGate_' baseName '.mat'],'file')==2 %for backwards compatibility with Bre gate files
					gateFileName = ['JCGate_' baseName '.mat'];
				end
				
			else
				gateFileName=gateFile;
			end

			GF = load(gateFileName);

			[fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(sampleFile);
			ch_SSC_A = getChannel(fcshdr,'SSC-A'); %Get channel number for SSC-A
			ch_FSC_A = getChannel(fcshdr,'FSC-A'); %Get channel number for FSC-A
			ch_FSC_W = getChannel(fcshdr,'FSC-W'); %Get channel number for FSC-W
			ch_FSC_H = getChannel(fcshdr,'FSC-H'); %Get channel number for FSC-H
			ch_SSC_W = getChannel(fcshdr,'SSC-W'); %Get channel number for SSC-W
			ch_SSC_H = getChannel(fcshdr,'SSC-H'); %Get channel number for SSC-H
			ch_time = getChannel(fcshdr,'Time'); %Get channel number for time
			ch_blue = getChannel(fcshdr,'Blue'); %Get channel number for Pacific-Blue
			ch_cyan = getChannel(fcshdr,'Cyan'); %Get channel number for AmCyan
			ch_yellow = getChannel(fcshdr,'FIT'); %Get channel number for FITC
			ch_red = getChannel(fcshdr,'Red'); %Get channel number for Texas-Red


			[P1inds g1] = Gating.gatePolygon(fcsdat(:,ch_FSC_A),fcsdat(:,ch_SSC_A),'semilogy',GF.gate1);
			[P2inds g2] = Gating.gatePolygon(fcsdat(:,ch_FSC_W),fcsdat(:,ch_FSC_H),'linear',GF.gate2);
			[P3inds g3] = Gating.gatePolygon(fcsdat(:,ch_SSC_W),fcsdat(:,ch_SSC_H),'semilogy',GF.gate3);
			JCInds = P1inds & P2inds & P3inds;

			if isfield(GF,'gate4') && isfield(GF,'gate5') && isfield(GF,'gate6') % modified by JG. only does this if time is applied in the gate file
				[P4inds g4] = Gating.gatePolygon(fcsdat(:,ch_time),fcsdat(:,ch_blue),'semilogy',GF.gate4);
				[P5inds g5] = Gating.gatePolygon(fcsdat(:,ch_time),fcsdat(:,ch_yellow),'semilogy',GF.gate5);
				[P6inds g6] = Gating.gatePolygon(fcsdat(:,ch_time),fcsdat(:,ch_red),'semilogy',GF.gate6);
				JCInds = P1inds & P2inds & P3inds & P4inds & P5inds & P6inds;
			end


			fcsdatGated = fcsdat(JCInds,:);

		end
		
		
		function createJCGate(fileName,use_time)
			% creates a file with gates for each of the Just Cell gates:
			%       gate1: SSC-A vs. FSC-A
			%       gate2: FSC-H vs. FSC-W
			%       gate3: SSC-H vs. SSC-W
			%       gate4: Time vs. Blue
			%       gate5: Time vs Yellow
			%       gate6: Time vs Red

			%decide if you want to use time gating as well as FSC/SSC

			%modified by JG to include time gating if user specifies
			%and also to use logicals when gating rather than indeces

			if nargin<2
				use_time = 0;
			end

			baseName=fileName(1:end-4);
			gateFileName = [baseName]; %modified by JG. previously gateFileName=['JCGate_' baseName];
			alternate_file_name = ['JCGate_' baseName]; %added for backwards compatibility with Bre's previous gate files


			if ~( exist([gateFileName '.mat'],'file')==2 || exist([alternate_file_name '.mat'],'file')==2 ) %only make a new gate file if one doesn't exist yet
				
				[fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(fileName);
				
				%P1 - Create and apply gate 1
				ch_SSC_A=getChannel(fcshdr,'SSC-A');
				ch_FSC_A=getChannel(fcshdr,'FSC-A');
				
				datSSC_A=fcsdat(:,ch_SSC_A);
				datFSC_A=fcsdat(:,ch_FSC_A);
				
				figure
				xlabel('FSC-A')
				ylabel('SSC-A')
				hold on
				[P1_logicals, gate1] = Gating.gatePolygon(datFSC_A, datSSC_A, 'semilogy');
				
				
				%P2 - Create and apply gate 2
				ch_FSC_W = getChannel(fcshdr, 'FSC-W');
				ch_FSC_H = getChannel(fcshdr, 'FSC-H');
				
				datFSC_W = fcsdat(:, ch_FSC_W);
				datFSC_W(~P1_logicals) = NaN;
				
				datFSC_H = fcsdat(:, ch_FSC_H);
				datFSC_H(~P1_logicals) = NaN;
				
				
				figure
				xlabel('FSC-W')
				ylabel('FSC-H')
				hold on
				[P2_logicals, gate2] = Gating.gatePolygon(datFSC_W, datFSC_H, 'linear');
				
				
				%P3 - Create and apply gate 3
				ch_SSC_W=getChannel(fcshdr,'SSC-W');
				ch_SSC_H=getChannel(fcshdr,'SSC-H');
				
				datSSC_W=fcsdat(:,ch_SSC_W);
				datSSC_W(~(P1_logicals & P2_logicals))= NaN;
				
				datSSC_H=fcsdat(:,ch_SSC_H);
				datSSC_H(~(P1_logicals & P2_logicals))= NaN;
				
				figure
				xlabel('SSC-W')
				ylabel('SSC-H')
				hold on
				[P3_logicals gate3] = Gating.gatePolygon(datSSC_W,datSSC_H,'semilogy');
				
				if use_time ~=0 %only apply time gating if user indicates so
					
					%P4 - Create and apply gate 4
					ch_x=getChannel(fcshdr,'FIT');
					ch_y=getChannel(fcshdr,'Red');
					
					dat_time=fcsdat(:,ch_x);
					dat_time(~(P1_logicals & P2_logicals & P3_logicals))= NaN;
					
					dat_blue=fcsdat(:,ch_y);
					dat_blue(~(P1_logicals & P2_logicals & P3_logicals))= NaN;
					
					
					figure
					xlabel('Time')
					ylabel('Pacific Blue')
					hold on
					[P4_logicals gate4] = Gating.gatePolygon(dat_time,dat_blue,'loglog');
					
					
					%P5  - Create and apply gate 5
					ch_x=getChannel(fcshdr,'Time');
					ch_y=getChannel(fcshdr,'FIT');
					
					dat_time=fcsdat(:,ch_x);
					dat_time(~(P1_logicals & P2_logicals & P3_logicals & P4_logicals))= NaN;
					
					dat_yellow=fcsdat(:,ch_y);
					dat_yellow(~(P1_logicals & P2_logicals & P3_logicals & P4_logicals))= NaN;
					
					
					figure
					xlabel('Time')
					ylabel('Yellow')
					hold on
					[P5_logicals gate5] = Gating.gatePolygon(dat_time,dat_yellow,'semilogy');
					
					
					%P6 - Create and apply gate 6
					ch_x=getChannel(fcshdr,'Time');
					ch_y=getChannel(fcshdr,'Red');
					
					dat_time=fcsdat(:,ch_x);
					dat_time(~(P1_logicals & P2_logicals & P3_logicals & P4_logicals & P5_logicals))= NaN;
					
					dat_red=fcsdat(:,ch_y);
					dat_red(~(P1_logicals & P2_logicals & P3_logicals & P4_logicals & P5_logicals))= NaN;
					
					figure
					xlabel('Time')
					ylabel('Red')
					hold on
					[P6_logicals gate6] = Gating.gatePolygon(dat_time,dat_red,'semilogy');
				end
				
				
				if use_time == 0 %Save gates 1-3
					%JCInds=P1_logicals & P2_logicals & P3_logicals; %JG-Don't think this line does anything here (?)
					save(gateFileName,'gate1','gate2','gate3')
					
				elseif use_time == 1 %Save gates 1-6
					%JCInds=P1_logicals & P2_logicals & P3_logicals & P4_logicals & P5_logicals & P6_logicals; %JG-Don't think this line does anything here (?)
					save(gateFileName,'gate1','gate2','gate3','gate4','gate5','gate6')
				end
				
			else
				warnmsg=[gateFileName '.mat already exhists'];
				warning(warnmsg)
			end

			hold off

		end
		
	end
	
end