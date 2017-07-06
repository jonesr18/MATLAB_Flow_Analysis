classdef Transforms < handle
	% Compilation of data transforms from the Weiss Lab flow cytometry MATLAB repository
	%
	%	Methods are implemented as static functions, so they can be called directly
	%	or after creating a handle to this class.
	%
	%		Ex:
	%		tfData = Transforms.lin2logicle(data);
	%
	%		Ex2:
	%		tf = Transforms();
	%		tfData = tf.logicle2lin(data);
	%
    %   Functions:
    %
    %       out = lin2logicle(S)
    %       out = logicle2lin(X)
    %       data = fcs2MEFL(data, channelFits, dataType)
    %       [channelFits, beadDir] = calibrateMEFL(beadsFilename)
    %       data = fcs2ABC(data, channelFits, dataType)
    %       [channelFits, locs, peaks] = calibrateABC(beadsFilename, channel, largePeaks, locs, peaks)
    %
    % Written/Compiled by Ross Jones
    % Weiss Lab, MIT
    % Last updated 2016-05-27
    
	methods (Static)
		
		function out = lin2logicle(S)
			% LIN2LOGICLE(S) converts a vector of values from a linear scale to a logicle scale 
			%   OUT = LIN2LOGICLE(S) returns a vector of values corresponding to the values 
			%   in X transformed from a linear scale to a logicle one using
			%   the conversion function described by Parks, et al. "A New Logicle
			%   Display Method Avoids Deceptive Effects of Logarithmic Scaling for Low
			%   Signals and Compensated Data" Equation (5), where
			%     S - a vector of linear 'raw' values
			%     out - a vector of logicle converted values from S
			%
			%   Since Equation 5 cannot be explicitly solved for X, a spline is fitted to
			%   X and S data using logicle2lin.m.  This spline is saved in the same folder 
			%   and used to evaluated the values in S.
			%
			%   Example:
			%       out = Transforms.lin2logicle(linspace(0,1000));
			%
			%   Written by
			%   Breanna Stillo
			%   bstillo@mit.edu
			%   Last Updated: 2014-10-14
			
			X = linspace(0, 4.5, 1000); % 4.5 = lin2logicle(2^18)
			Y = Transforms.logicle2lin(X);
			p = spline(Y, X);

			out = ppval(p, S);

		end
		
		
		function out = logicle2lin(X)
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
			%       out = Transforms.logicle2lin(linspace(0,4));
			%
			%   Written by
			%   Breanna Stillo
			%   bstillo@mit.edu
			%   Last Updated: 2014-10-14
			%
			%   Edit 2015-02-09 (Ross) - Added comments and intuitive variable names
			%                            Major speed up by using algebraic conversion 
            %                            instead of symbolic toolbox
            %   Edit 2016-03-28 (Ross) - Made display breadth smaller w/ Bre's adjustment of 
            %                            DATA_MIN, which increases the relative size of the area
            %                            between 10^-2 and 10^2

			% Set values
			DATA_MAX = 2^18;         % (T) The maximum measurable data value ( = 26214 )
			DISPLAY_BREADTH = 4.5;   % (M) Breadth of the display in decades ( = lin2logicle(2^18) )
			DATA_MIN = -150;         % (r) Negative range reference value    ( lin2logicle(0) = -171 )

			% This gives a range for linearization around zero (W).
			linRange = (DISPLAY_BREADTH - log10(DATA_MAX / abs(DATA_MIN))) / 2;

			% p and W (linRange) are considered one parameter, p is introduced by the 
			% authors for compactness. Here we find p that solves their equivalence:
			%   w = 2p ln(p)/(p + 1)
			% Solved by WolframAlpha:
			P = linRange / (2 * lambertw(0.5 * exp(-linRange / 2) * linRange));

			% Find where the given data is out of the linear range.
			logX = (X >= linRange);

			% Compute and return the final linearized vector 
			out = toLin(X) .* logX - toLin(2 * linRange - X) .* (1 - logX);
			
			
			% --- Helper Functions --- %
			
			
			function converted = toLin(convert)
				% Does a final conversion
				
				converted = DATA_MAX .* 10^-(DISPLAY_BREADTH - linRange) .* ...
					(10.^(convert - linRange) - P^2 .* 10.^(-(convert - linRange) ./ P) + P^2 - 1);
			end
        end
		
        
		function data = fcs2MEFL(data, channelFits, dataType)
			%out = fcs2MEFL(data, channelFits)
			% converts fcs data to MEFL
            %
            %   Inputs
            %
			%       data            A standard data struct. The channels must match those in 
            %                       channelFits. A new field called mefl is added for each 
            %                       converted (non-zero slope) channel containing the MEFL data.
            %
            %       channelFits     (Output from calibrateMEFL())
            %                       An Nx2 table of slopes/intercepts for linear conversion of 
            %                       N channels to MEFL equivalents. Rows are labeled with 
            %                       channel names and columns are labeled with slope / zero
            %
            %       dataType        The population to convert to MEFLs (eg 'scComp' or 'tf')
            %                        - must be in data!
			%
			%   Written by
			%   Breanna Stillo
			%   bstillo@mit.edu
			%
            %   Update 2016-04-23 by Ross Jones
            %       Uses standard struct interface and uses channelFits directly rather than
            %       loading the saveFile
			
            % Process inputs
			channels = checkInputs(data, channelFits, dataType);
            
			r = 2^18; %resolution
			n = 4.5; %log decades
            
			for i = 1:numel(data) %#ok<ALIGN>
                
                for ch = 1:numel(channels)
                
                    channel = channels{ch};
                    chData = data(i).(channel).(dataType);
                    
                    % Remove negative values
                    if any(chData < 0)
                        warning('Negative values encountered, removing')
                        chData = chData(chData >= 0);
                    end
                    
                    % Do conversion
                    m = channelFits(channel, :).slope;
                    b = channelFits(channel, :).zero;
                    data(i).(channel).mefl = 10^b .* chData.^(m * r / n); 
                    
                end
            end
			
			
			% --- Helper Functions --- %
			
            
			function channels = checkInputs(data, channelFits, dataType)
				
				validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                validateattributes(channelFits, {'table'}, {}, mfilename, 'channelFits', 2);
                
                tableChannels = channelFits.Properties.RowNames;
                dataChannels = data(1).chanNames;
                if ~isempty(setdiff(tableChannels, dataChannels))
                    error('Channels in channelFits do not match channels in data!')
                end
                
                channels = tableChannels;
                validatestring(dataType, fieldnames(data(1).(channels{1})), mfilename, 'dataType', 3);
			end

        end
		
        
        function [channelFits, beadDir] = calibrateMEFL(beadsFilename)
            % calibrateMEFL creates MEFL fits for an .fcs file of a particular bead
            % sample chosen by the user. Note this is for rainbow calibration beads, not 
            % other antibody-tagged beads (see calibrateABC()). 
            %
            % A subfolder is created with the name of the bead file where associated calibration 
            % plots are stored and the calibration table corresponding to MEFL fits.
            %
            %   Inputs: 
            %   
            %       beadFilename        An .fcs filename which contains rainbow bead data
            %                           If a cell array of strings, the data from each file 
            %                           will be concatenated together.
            %
            %   Outputs:
            %
            %       channelFits         An Nx2 table of slopes/intercepts for linear conversion of 
            %                           N channels to MEFL equivalents. Rows are labeled with 
            %                           channel names and columns are labeled with slope / zero
            %       
            %       beadDir             The name of the directory where bead fits are stored
            %
            %   Written by
            %   Breanna Stillo
            %   bstillo@mit.edu
            %   Last Updated: 2016-04-23 by Ross Jones

            % Iterate over all given filenames and concatenate data together. 
            if ~iscell(beadsFilename)
                beadsFilename = {beadsFilename};
            end
            data = FlowAnalysis.openFiles(beadsFilename{:});

            beadDir = strcat('calibration_', beadsFilename{1}(1:end-4));

            if exist(beadDir, 'dir')
                error('Directory %s already exists!\n', beadDir);
            else

                mkdir(beadDir)
                %channels not to convert to MEFL
                ignore={'Time','Event','SSC','FSC'};
                
                % Initialize as an array, convert to a table later
                channels = data(1).chanNames;
                channelFits = zeros(numel(channels), 2);
               
                % Iterate over channels and update as needed
                for chID = 1:numel(channels)
                    
                    channel = channels{chID};
                    fprintf(1, 'Processing %s\n', channel);

                    % Check if name is in the 'ignore' list
                    skip = false;
                    for ign=1:length(ignore)  
                        if strfind(channel, ignore{ign})
                            skip = true;
                            break
                        end
                    end
                    
                    if skip
                        channelFits(chID, :) = [1 0];
                        fprintf('skipped\n')
                    else
                        
                        % Extract channel data
                        channelData = [];
                        for i = 1:length(data)
                            channelData = [channelData; data(i).(channel).raw]; %#ok<AGROW>
                        end
                        numBeads = length(channelData);
                        
                        figure('Position', [230 50 350 600])
                        ax1 = subplot(2, 1, 1);
                        hold(ax1, 'on');
                        numBins = round(numBeads / 100);
                        [nelements, centers] = Plotting.biexhist(channelData, numBins);

                        % Find highest peaks
                        peaks = findpeaks(nelements);
                        maxpeak = max(peaks);
                        [peaks, locs] = findpeaks(nelements, 'MinPeakHeight', maxpeak / 10);
                        plot(ax1, centers(locs), peaks, '*r')
                        title(strrep(channel, '_', '-'))
                        ylabel('COUNT')

                        % Assume we have the P highest peaks
                        P = min(length(peaks), 7);
                        MEFL = [792 2079 6588 16471 47497 137049 271647];
                        MEFL = MEFL(end - P + 1 : end);
                        peakIntensity = Transforms.logicle2lin(centers(locs));
                        peakIntensity = peakIntensity(end - P + 1 : end);

                        fit = polyfit(log10(peakIntensity), log10(MEFL), 1);

                        ax2 = subplot(2, 1, 2);
                        plot(ax2, ...
                            log10(peakIntensity), MEFL, 'ob', ...
                            log10(peakIntensity), 10.^polyval(fit, log10(peakIntensity)), 'r-')
                        ax2.XScale = 'log';
                        ax2.YScale = 'log';
                        ax2.XLim = [1e1, 1e5];
                        ax2.YLim = [1e2, 1e7];
                        ax2.FontSize = 14;
                        
                        
                        yresid = MEFL - 10.^polyval(fit, log10(peakIntensity));
                        SSresid = sum(yresid.^2);
                        SStotal = (length(MEFL) - 1) * var(MEFL);
                        rsq = 1 - SSresid / SStotal;
                        title(sprintf('R^2 = %.3f', rsq))
                        xlabel('Fluorescence')
                        ylabel('MEFL')
                        hold off
                        
                        % Assign fit data
                        channelFits(chID, :) = fit;

                        fprintf(1, 'Finished fitting %s. # Peaks = %d\n', channel, length(peaks));

                        saveas(gcf, [beadDir '\' channel '_MEFL_Fit.fig'])
                    end
                end
                
                % Convert fits to a table
                channelFits = array2table(channelFits);
                channelFits.Properties.VariableNames = {'slope', 'zero'};
                channelFits.Properties.RowNames = channels;
                
                % Save channelFits
                save([beadDir '\' beadDir '_MEFL_Fit.mat'], 'channelFits')
            end
        end
        
        
        function data = fcs2ABC(data, channelFits, dataType)
			%out = fcs2MEFL(data, channelFits)
			% converts fcs data to MEFL
            %
            %   Inputs
            %
			%       data            A standard data struct. The channels must match those in 
            %                       channelFits. A new field called mefl is added for each 
            %                       converted (non-zero slope) channel containing the MEFL data.
            %
            %       channelFits     (Output from calibrateABC())
            %                       An Nx2 table of slopes/intercepts for linear conversion of 
            %                       N channels to MEFL equivalents. Rows are labeled with 
            %                       channel names and columns are labeled with slope / zero
            %
            %       dataType        The population to convert to MEFLs (eg 'scComp' or 'tf')
            %                        - must be in data!
			%
			%   Written by
			%   Ross Jones
			%   jonesr18@mit.edu
			%   2016-04-23
            
            % Process inputs
			channels = checkInputs(data, channelFits, dataType);
                        
			for i = 1:numel(data) %#ok<ALIGN>
                
                for ch = 1:numel(channels)
                
                    channel = channels{ch};
                    chData = data(i).(channel).(dataType);
                    
%                     % Remove negative values -- might be better to allow isolation post-proc
%                     if any(chData < 0)
%                         warning('Negative values encountered, removing')
%                         chData = chData(chData >= 0);
%                     end
                    
                    % Do conversion
                    m = channelFits(channel, :).slope;
                    b = channelFits(channel, :).zero;
                    data(i).(channel).abc = 10^b .* chData.^m;                    
                end
            end
			
			
			% --- Helper Functions --- %
			
            
			function channels = checkInputs(data, channelFits, dataType)
				
				validateattributes(data, {'struct'}, {}, mfilename, 'data', 1);
                validateattributes(channelFits, {'table'}, {}, mfilename, 'channelFits', 2);
                
                tableChannels = channelFits.Properties.RowNames;
                dataChannels = data(1).chanNames;
                if ~isempty(setdiff(tableChannels, dataChannels))
                    error('Channels in channelFits do not match channels in data!')
                end
                
                channels = tableChannels;
                validatestring(dataType, fieldnames(data(1).(channels{1})), mfilename, 'dataType', 3);
			end

        end
        
        
        
        function channelFits = calibrateABCfromMedians(beadsFilenames, channel, hasBlank)
            % calibrateABCfromMedians creates antibody binding capacity (ABC) fits for an 
            % .fcs file of a particular bead sample supplied by the user using medians of 
            % each file. Note this is for Quantum Simply Cellular (R) beads, not rainbow 
            % beads (see calibrateMEFL()). 
            %
            % Calibration plots and the table with ABC fits are saved in pwd().
            %
            %   Inputs: 
            %   
            %       beadFilename        A cell array of .fcs filenames which contain Ab-stained 
            %                           bead data. Can be given in any order (median is sorted).
            %
            %       channel             The channel name to calibrate bead data
            %
            %       hasBlank            (Optional) Boolean to indicate if the samples include the
            %                           blank population (ignore lowest peak if so). Default = false
            %
            %   Outputs:
            %
            %       channelFit          An Nx2 table of slopes/intercepts for linear conversion of 
            %                           N channels to MEFL equivalents. Rows are labeled with 
            %                           channel names and columns are labeled with slope / zero
            %       
            %       medians             A 1xN array of medians values in order the files are given
            %
            % Written by Ross Jones
            % 2016-06-01
            % MIT Weiss Lab
            % jonesr18@mit.edu
            
            % Import data
            beadData = FlowAnalysis.openFiles(beadsFilenames{:});
            
            % Initialize as an array, convert to a table later
            channelFits = zeros(1, 2);
            
            % Calculate medians
            beadMedians = zeros(1, numel(beadData));
            for i = 1:numel(beadData)
                % Check if unimodal (if not, there are many unstained beads)
%                 sortedBeads = sort(beadData(i).(channel).raw);
%                 logSortedBeads = Transforms.lin2logicle(sortedBeads);
%                 [dip, p_value] = hartigansDipSignifTest(logSortedBeads, 5000);
%                 if (p_value < 0.05)
%                     warning('Significant non-unimodal behavior detected: Dip = %.2f, P = %2f', ...
%                         dip, p_value);
%                     
%                     % If multi-modal, compute mean of higher population (stained beads) with
%                     % gaussian mixture model
%                     figure()
%                     Plotting.biexhist(sortedBeads);
%                     gmm = fitgmdist(logSortedBeads, 2, 'Start', 'plus', 'Replicates', 5);
%                     beadMedians(i) = Transforms.logicle2lin(max(gmm.mu));
%                 else
                    beadMedians(i) = median(beadData(i).(channel).raw);
%                 end
            end
            beadMediansSorted = sort(Transforms.lin2logicle(beadMedians));
            
            % Extract channel data and plot as histogram
            figPeaks = figure();
            ax1 = subplot(1, 2, 1);
            hold(ax1, 'on');
            for i = 1:length(beadData)
                channelData = beadData(i).(channel).raw; 
                numBeads = length(channelData);
                numBins = numBeads / 100;
                [binCounts, binCenters] = Plotting.biexhist(channelData, numBins);
                
                % Match closest bin to medians
                [~, binIdx] = min(abs(binCenters - Transforms.lin2logicle(beadMedians(i))));
                peak = binCounts(binIdx);
                plot(ax1, Transforms.lin2logicle(beadMedians(i)), peak, '*r')
            end         
            title('Bead Peak Medians')
            ylabel('COUNT')
            xlabel(strrep(channel, '_', '-'))
            ax1.FontSize = 14;
            hold(ax1, 'off')
            
            % Assume we have the P highest peaks and an unstained population
            ABC = [5899, 68458, 281582, 621651];
            if (exist('hasBlank', 'var') && hasBlank)
                P = min(length(beadMediansSorted) - 1, numel(ABC));
            else
                P = min(length(beadMediansSorted), numel(ABC));
            end
            ABC = ABC(end - P + 1 : end);
            peakIntensity = Transforms.logicle2lin(beadMediansSorted);
            peakIntensity = peakIntensity(end - P + 1 : end);
            
            % Calculate fit
            fit = polyfit(log10(peakIntensity), log10(ABC), 1);
            
            ax2 = subplot(1, 2, 2);
            plot(ax2, ...
                peakIntensity, ABC, 'ob', ...
                peakIntensity, 10.^polyval(fit, log10(peakIntensity)), 'r-')
            ax2.XScale = 'log';
            ax2.YScale = 'log';
            ax2.XLim = [1e-1, 1e5];
            ax2.YLim = [1e3, 1e7];
            ax2.FontSize = 14;
            
            % Calculate R^2 value
            yresid = ABC - 10.^polyval(fit, log10(peakIntensity));
            SSresid = sum(yresid.^2);
            SStotal = (length(ABC) - 1) * var(ABC);
            rsq = 1 - SSresid / SStotal;
            title(sprintf('R^2 = %.3f', rsq))
            xlabel('Fluorescence')
            ylabel('ABC')
            hold(ax2, 'off')

            % Assign fit data
            channelFits(1, :) = fit;

            fprintf(1, 'Finished fitting %s. # Peaks = %d\n', channel, length(peaks));
            savefig(figPeaks, [channel '_ABC_Fit.fig'])

            % Convert fits to a table
            channelFits = array2table(channelFits);
            channelFits.Properties.VariableNames = {'slope', 'zero'};
            channelFits.Properties.RowNames = {channel};

            % Save channelFits
            save([channel, '_ABC_Fit.mat'], 'channelFits')
        end
        
        
        
        function [channelFits, locs, peaks] = calibrateABC(beadsFilename, channel, largePeaks, hasBlank, locs, peaks)
            % calibrateABC creates antibody binding capacity (ABC) fits for an .fcs file 
            % of a particular bead sample supplid by the user. Note this is for Quantum 
            % Simply Cellular (R) beads, not rainbow beads (see calibrateMEFL()). 
            %
            % Calibration plots and the table with ABC fits are saved in pwd().
            %
            %   Inputs: 
            %   
            %       beadFilename        An .fcs filename which contains Ab-stained bead data
            %                           If a cell array of strings, the data from each file 
            %                           will be concatenated together
            %
            %       channel             The channel name to calibrate bead data
            %
            %       largePeaks          (Optional) Boolean to indicate whether to rule out peaks
            %                           smaller than 1/10 max(peaks) on the first peak finding 
            %                           run. Default = true
            %
            %       hasBlank            (Optional) Boolean to indicate if the samples include the
            %                           blank population (ignore lowest peak if so). Default = false
            %
            %       locs                (Optional) 1xN array of peak locations in ascending order
            %
            %       peaks                 (Optional) 1xN array of peak heights
            %
            %   Outputs:
            %
            %       channelFit          An Nx2 table of slopes/intercepts for linear conversion of 
            %                           N channels to MEFL equivalents. Rows are labeled with 
            %                           channel names and columns are labeled with slope / zero
            %
            %       locs                A 1xN array of peak locations in ascending order
            %
            %       peaks                 A 1xN array of peak heights
            %
            % Written by Ross Jones
            % 2016-04-23
            % MIT Weiss Lab
            % jonesr18@mit.edu
            % 
            % Updated 2016-06-01
            
            % Iterate over all given filenames and concatenate data together. 
%             assert(iscell(beadsFilenames), 'beadsFilename must be a cell array of strings!');
            if ~iscell(beadsFilename)
                beadsFilename = {beadsFilename};
            end
            beadData = FlowAnalysis.openFiles(beadsFilename{:});
            if ~exist('largePeaks', 'var')
                largePeaks = true;
            end
            
            % Initialize as an array, convert to a table later
            channelFits = zeros(1, 2);

            % Extract channel data
            channelData = [];
            for i = 1:length(beadData)
                channelData = [channelData; beadData(i).(channel).raw]; %#ok<AGROW>
            end
            numBeads = length(channelData);
            
            % Plot figure to show fitting
            figPeaks = figure('Position', [500 500 700 350]);
            subplot(1, 2, 1)
            numBins = round(numBeads / 100);
            [counts, centers] = Plotting.biexhist(channelData, numBins);
            hold on

            % Find highest peaks
            if ~(exist('locs', 'var') && exist('peaks', 'var'))
                [peaks, locs] = findpeaks(counts);
                if (largePeaks)
                    maxpeak = max(peaks);
                    [peaks, locs] = findpeaks(counts, 'MinPeakHeight', maxpeak / 10);
                end
            end
            
            % Plot peaks on data
            subplot(1, 2, 1)
            hold on
            plot(centers(locs), peaks, '*r')
            title(strrep(channel, '_', '-'))
            ylabel('COUNT')
            
            % Assume we have the P highest peaks and an unstained population
            ABC = [68458, 281582, 621651];
%             ABC = [5899, 68458, 281582, 621651];

            if (exist('hasBlank', 'var') && hasBlank)
                P = min(length(locs) - 1, numel(ABC));
            else
                P = min(length(locs), numel(ABC));
            end
            ABC = ABC(end-P+1:end);
            peakIntensity = Transforms.logicle2lin(centers(locs));
            peakIntensity = peakIntensity(end-P+1:end);
            
            % Calculate fit
            r = 2^18;   % resolution
            n = 5;      % log decades
            peakRelChannel = r / n * log10(peakIntensity);
            fit = polyfit(peakRelChannel, log10(ABC), 1);
            
            subplot(1,2,2)
            semilogy(peakRelChannel, ABC, 'ob')
            hold on
            semilogy(peakRelChannel, 10.^polyval(fit, peakRelChannel), 'r-')
            
            % Calculate R^2 value
            yresid = ABC - 10.^polyval(fit, peakRelChannel);
            SSresid = sum(yresid.^2);
            SStotal = (length(ABC) - 1) * var(ABC);
            rsq = 1 - SSresid / SStotal;
            title(sprintf('R^2 = %.3f', rsq))
            axis([0 2^18 1 1e6])
            xlabel('REL CH #')
            ylabel('ABC')
            hold off

            % Assign fit data
            channelFits(1, :) = fit;

            fprintf(1, 'Finished fitting %s. # Peaks = %d\n', channel, length(peaks));
            savefig(figPeaks, [channel '_ABC_Fit.fig'])

            % Convert fits to a table
            channelFits = array2table(channelFits);
            channelFits.Properties.VariableNames = {'slope', 'zero'};
            channelFits.Properties.RowNames = {channel};

            % Save channelFits
            save([channel, '_ABC_Fit.mat'], 'channelFits')
        end
        
    end
    
end