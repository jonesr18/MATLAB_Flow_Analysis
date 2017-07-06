function [xBleed,yComp]=compensateSC_1028(JCdat, SCdat, sampledat, plotsOn)
%[xBleed,yComp]=compensate(JCFile, singleColorFile, correctFile, bC, cC,plotsOn)
%single-color compensation
    %JCFile = Just Cell File
    %singleColorFile = FCS file for single color control
    %correctFile = FCS file that needs to be corrected
    %bC = channel number of the color that is bleeding and needs to be
    %     corrected for in other channels
    %cC = channel number of the color that needs to be corrected
    %plotsOn = optional. Turns plots on for visualization
    
%   Example:
%       [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs('sample.fcs');
%       GreenChannel = getChannel(fcshdr,'FIT');
%       RedChannel = getChannel(fcshdr,'Red');
%       GreenData = fcsdat(:,GreenChannel);
%       RedData = fcsdat(:,RedChannel);
%       figure
%       biexplot(redData,greenData)
%       [newRed,newGreen]=compensateSC('JCsample.fcs','RedSample.fcs','sample.fcs',RedChannel,GreenChannel);
%       figure
%       biexplot(newred,newgreen)
%
%   Written by
%   Breanna Stillo
%   bstillo@mit.edu
%   Last Updated: 2014-10-21;

%     [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(JCFile);
    cVals = JCdat(:,2);
    
    cJC = median(cVals);

    
%     [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(singleColorFile);

    bVals = SCdat(:,1);
    cVals = SCdat(:,2);
    
    XI=[0 10.^(linspace(2,log10(max(bVals)),20))];
    YI =lsq_lut_piecewise(bVals,cVals,XI);
    
    if exist('plotsOn','var')
        if plotsOn==1
            figure
            biexplot(bVals,cVals)
            hold on
            biexplot(XI,YI,'r-')
        end
    end
    
%     [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(correctFile);

    bVals = sampledat(:,1);
    cVals = sampledat(:,2);
    
    yerror=compEval(bVals);
    cReal=cVals-yerror+cJC;
    
    xBleed=bVals;
    yComp=cReal;
    
    function y=compEval(x)
                
        yBuff=zeros(length(x),1);
        for w=1:length(x)
            eqn=find(XI>x(w),1,'first');
            if isempty(eqn)
                eqn=length(XI);
            end
            if eqn==1
                eqn=2;
            end

            a=(YI(eqn)-YI(eqn-1))./(XI(eqn)-XI(eqn-1));
            b=YI(eqn)-a*XI(eqn);
            yBuff(w)=a*x(w)+b;
        end
        
        y=yBuff;
            
    end



end