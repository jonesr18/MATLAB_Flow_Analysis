function [xBleed,yComp]=bleedSC(JCFile, singleColorFile, correctFile, bC, cC,plotsOn)
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
%   Edit 11/20/2015 Ross Jones:
%   replace linear piecewise fitting with 'interp1' for smoothing and speed
%
%   Edit 12/10/2015 Bre Stillo:
%   Add buffers to ends of XI, YI to accomodate new logicleTransform
%   function
%
%   Written by
%   Breanna Stillo
%   bstillo@mit.edu
%   Last Updated: 2014-10-21;

%     [fcsdatJC, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(JCFile);
    fcsdatJC = applyJCGate(JCFile, JCFile);
    cValsJC = lin2logicle(fcsdatJC(:,cC));
    bValsJC = lin2logicle(fcsdatJC(:,bC));
%     cValsJC = fcsdatJC(:,cC);
%     bValsJC = fcsdatJC(:,bC);
    
    clow = prctile(cValsJC,1);
    chigh = prctile(cValsJC,99);
    blow = prctile(bValsJC,1);
    bhigh = prctile(bValsJC,99);
%     cJC = mean(cValsJC);
%     bJC = mean(bValsJC);
%     cJC = median(cValsJC);
%     bJC = median(bValsJC);

    
%     [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(singleColorFile);
    fcsdat = applyJCGate(singleColorFile,JCFile);

    bVals = lin2logicle(fcsdat(:,bC));
    cVals = lin2logicle(fcsdat(:,cC));
    
    cJC = mean(cVals(cVals>clow & cVals<chigh));
    bJC = mean(bVals(bVals>blow & bVals<bhigh));
%     bVals = fcsdat(:,bC);
%     cVals = fcsdat(:,cC);
    
%     M = sortrows([bVals cVals],1);
%     bVals = M(:,1);
%     cVals = M(:,2);
%     
%     cVals = smooth(cVals,100);
%     cVals = smooth(cVals);

    numPosPoints = sum(bVals>bhigh);
    numBins = round(numPosPoints/600);
    XI = [linspace(prctile(bValsJC,99),prctile(bVals,98),numBins)];
    YI = lsq_lut_piecewise(bVals,cVals,XI);
    YI(1) = cJC;
%     XI = [lin2logicle(-1000) bJC XI XI(end)+XI(end)-XI(end-1)];
    XI = [lin2logicle(-1000) bJC XI max(bVals) max(bVals)+(max(bVals)-XI(end))];
    YI = [YI(1); cJC; YI; YI(end)+ (YI(end)-YI(end-1))/(XI(end-2)-XI(end-3))*(max(bVals)-XI(end-2));YI(end)+ 2*(YI(end)-YI(end-1))/(XI(end-2)-XI(end-3))*(max(bVals)-XI(end-2))];
    
    if exist('plotsOn','var')
        if plotsOn==1
            figure
            plot(bVals,cVals,'.')
            hold on
            plot(XI,YI,'ro-','LineWidth',2)
            plot(bJC,cJC,'b*')
        end
    end
    
%     YI = logicle2lin(YI);
    
%     [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(correctFile);
    fcsdat = applyJCGate(correctFile,JCFile);

    bVals = lin2logicle(fcsdat(:,bC));
    cVals = lin2logicle(fcsdat(:,cC));
%     cVals = fcsdat(:,cC);
    
    yerror=interp1(XI,YI,bVals);
    cReal=cVals-yerror+cJC;
%     cReal=cVals-yerror+logicle2lin(cJC);
    
    xBleed=logicle2lin(bVals);
    yComp=logicle2lin(cReal);
%     yComp = cReal;
    
    
    
    
%     function y=compEval(x)
%                 
%         yBuff=zeros(length(x),1);
%         for w=1:length(x)
%             eqn=find(XI>x(w),1,'first');
%             if isempty(eqn)
%                 eqn=length(XI);
%             end
%             if eqn==1
%                 eqn=2;
%             end
% 
%             a=(YI(eqn)-YI(eqn-1))./(XI(eqn)-XI(eqn-1));
%             b=YI(eqn)-a*XI(eqn);
%             yBuff(w)=a*x(w)+b;
%         end
%         
%         y=yBuff;
%             
%     end



end