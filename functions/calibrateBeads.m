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

    % channel=getChannel(fcshdr,'FIT');
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
            numcells=length(channelData);

            figure('Position', [230 50 350 600])
            subplot(2,1,1)
            numBins=round(numcells/100);
            [nelements,centers]=biexhist(channelData,numBins);
            hold on

            %find highest peaks
            pks=findpeaks(nelements);
            maxpeak=max(pks);
            [pks,locs]=findpeaks(nelements,'MinPeakHeight',maxpeak/10);
            subplot(2,1,1)
            hold on
            plot(centers(locs),pks,'*r')
            title(ChanName)
            ylabel('COUNT')

            %assume we have the P highest peaks
            P=min(length(pks),7);
            MEFL=[792 2079 6588 16471 47497 137049 271647];
            MEFL=MEFL(end-P+1:end);
            peakIntensity=logicle2lin(centers(locs));
            peakIntensity=peakIntensity(end-P+1:end);

            r=2^18; %resolution
            n=5; %log decades
            peakRelChannel= r/n*log10(peakIntensity);

            fit=polyfit(peakRelChannel,log10(MEFL),1);
            
            subplot(2,1,2)
            semilogy(peakRelChannel,MEFL,'ob')
            hold on
            semilogy(peakRelChannel,10.^polyval(fit,peakRelChannel),'r-')
            
            yresid=MEFL-10.^polyval(fit,peakRelChannel);
            SSresid = sum(yresid.^2);
            SStotal = (length(MEFL)-1) * var(MEFL);
            rsq = 1 - SSresid/SStotal;
            t=['Rsqd = ' num2str(rsq)];
            title(t)
            axis([0 2^18 1 1e6])
            xlabel('REL CH #')
            ylabel('MEFL')
            hold off

            ChannelFits(channel,:) = fit;
            
            numPeaks=length(pks);
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

