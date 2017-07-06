% calibrateBeads creates MEFL fits for an fcs file of a particular bead
% sample chosen by the user.  A folder is created with the name of that
% bead file where associated calibration plots are stored and the
% calibration matrix corresponding to MEFL fits
%
%   Written by
%   Breanna Stillo
%   bstillo@mit.edu
%   Last Updated: 2014-10-22;

close all; clear all; clc;



beads_FCS = uigetfile('*.fcs','Select the Beads File for this experiement');

[fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(beads_FCS);

BeadFitName=['calibration_' beads_FCS(1:end-4)];

if exist([BeadFitName])~=7

    mkdir(BeadFitName)
    %channels not to convert to MEFL
    ignore={'Time','Event','SSC','FSC'};

    ChannelFits=zeros(length(fcsdat(1,:)),2);

    for channel=1:length(fcsdat(1,:))

        ChanName=fcshdr.par(channel).name;
        fprintf([ChanName '.......'])

        %check if name is in the 'ignore' list
        ok=1;
        for ign=1:length(ignore)  
            if strfind(ChanName,ignore{ign})
                ok=0;
            end
        end

        if ok        

            channelData=fcsdat(:,channel);
            channelData = channelData(channelData>0);
            rangeData = max(real(log10(channelData)))-min(real(log10(channelData)));
            numBins = round(rangeData/.02);
            
            [nelements, edges] = histcounts(real(log10(channelData)),numBins);
            centers = mean([edges(1:end-1); edges(2:end)],1);
            centers = 10.^centers;
            nelements = smooth(nelements,3);
            figure('Position', [230 50 350 600])
            voltage = fcshdr.par(channel).voltage;
            subplot(2,1,1)
            plot(centers,nelements,'k-','LineWidth',2)
            hold on
            area(centers,nelements,'FaceColor',[0.5 0.5 0.5])
            set(gca,'xscale','log','xlim',[1 1e6])
            

%             [nelements,centers]=biexhist(channelData,numBins);
%             hold on

            %find highest peaks
            nelements = smooth(nelements);
            pks=findpeaks(nelements);
            maxpeak=max(pks);
            [pks,locs]=findpeaks(nelements,'MinPeakHeight',maxpeak/10);%,'MinPeakProminence',maxpeak/6);
            
            NP = input(['\nHow many peaks for ' ChanName '? (Max of 7)\n']);
            if NP>7
                NP = input(['\nHow many peaks for ' ChanName '? (Max of 7!!!!!!!)\n']);
            elseif NP>locs
                NP = input(['\nNot enough peaks found, select lower number of peaks\nHow many peaks for ' ChanName '? (Max of 7)\n']);
            end
            
            toplocs = locs(end-(NP-1):end);
            
            subplot(2,1,1)
            hold on
%             plot(centers(locs),pks,'*r')
            for i = 1:length(toplocs)
                plot([centers(toplocs(i)), centers(toplocs(i))],get(gca,'ylim'),'-r','linewidth',2)
            end
            ylabel('COUNT')
            title([ChanName ' : Voltage = ' num2str(voltage)])

            %assume we have the P highest peaks
%             P=min(length(toplocs),7);
            MEFL=[792 2079 6588 16471 47497 137049 271647];
            MEFL=MEFL(end-NP+1:end);
            peakIntensity=centers(toplocs);
            peakIntensity=peakIntensity(end-NP+1:end);

%             r=2^18; %resolution
%             n=5; %log decades
%             peakRelChannel= r/n*log10(peakIntensity);

            fit=polyfit(log10(peakIntensity),log10(MEFL),1);
            
            subplot(2,1,2)
            loglog(peakIntensity,MEFL,'ob')
            hold on
            loglog(peakIntensity,10.^polyval(fit,log10(peakIntensity)),'r-')
            
            yresid=MEFL-10.^polyval(fit,log10(peakIntensity));
            SSresid = sum(yresid.^2);
            SStotal = (length(MEFL)-1) * var(MEFL);
            rsq = 1 - SSresid/SStotal;
            t=['Rsqd = ' num2str(rsq)];
            title(t)
            axis([1 1e6 1 1e6])
            xlabel(ChanName)
            ylabel('MEFL')
            hold off

            ChannelFits(channel,:) = fit;
            
            numPeaks=length(toplocs);
            printstr=['done. #peaks = ' num2str(numPeaks) '\n'];
            fprintf(printstr)
            
            saveas(gcf, [pwd '/' BeadFitName '/' ChanName '.fig'])
        else
            ChannelFits(channel,:)=[1 0];
            fprintf('skipped\n')
        end
    end
    
    eval([BeadFitName  '=ChannelFits;']);
    save([BeadFitName '/' BeadFitName],BeadFitName)
else
    errstr=[BeadFitName ' already exists\n'];
    error(errstr)
end

