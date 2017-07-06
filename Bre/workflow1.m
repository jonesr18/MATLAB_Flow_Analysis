function data = workflow1(sample_file, NT_file, SC_file, SC_ch)

save_name = ['workflow1/' sample_file(1:end-4) '.mat'];
plotting = 0;
if exist(save_name)
    return
else
    %% autofluorecence correction
    % calcualte autofluorescence value for each channel
    datNT = applyJCGate(NT_file, NT_file);
    autofluorescence = mean(datNT);
    
    % subtract autofluorescence
    datNT = datNT-repmat(autofluorescence,size(datNT,1),1);
    
    datS = applyJCGate(sample_file, NT_file);
    datS = datS-repmat(autofluorescence,size(datS,1),1);

    datSC = applyJCGate(SC_file,NT_file);
    datSC = datSC-repmat(autofluorescence,size(datSC,1),1);
    
    %% compensate1
%     [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs('s_R1.fcs');
    
    datNT_b = lin2logicle(datNT(:,SC_ch));
    datSC_b = lin2logicle(datSC(:,SC_ch));
    datS_b = lin2logicle(datS(:,SC_ch));
    
    datS_new = datS;
    
    % replace data in only color channels
    for ch=7:12
        if ch ~= SC_ch

            bhigh = prctile(datNT_b,99);

            datSC_c = lin2logicle(datSC(:,ch));
            datS_c = lin2logicle(datS(:,ch));

            numPoints = sum(datSC_b>bhigh);
            numBins = max(round(numPoints/800),5);
            XI = [linspace(bhigh,prctile(datSC_b,99),numBins)];
            YI = lsq_lut_piecewise(datSC_b,datSC_c,XI);
            XI = XI(2:end);
            YI = YI(2:end);
            XI = [lin2logicle(-1000) lin2logicle(0) XI max(datSC_b) max(datSC_b)+(max(datSC_b)-XI(end))];
            YI = [lin2logicle(0); lin2logicle(0); YI; YI(end)+ (YI(end)-YI(end-1))/(XI(end-2)-XI(end-3))*(max(datSC_b)-XI(end-2));YI(end)+ 2*(YI(end)-YI(end-1))/(XI(end-2)-XI(end-3))*(max(datSC_b)-XI(end-2))];
            YI_upper = repmat(4.5,length(YI),1);

            yerror = interp1(XI,YI,datS_b);
            top = 4.5;
            datS_c_new = top-(top-datS_c)./(top-yerror).*(top-lin2logicle(0));
            
            datS_new(:,ch) = logicle2lin(datS_c_new);
            
            if plotting
                figure
                subplot(2,2,1)
                hold on
                plot(datSC_b,datSC_c,'.')
                plot(XI,YI,'ro-','LineWidth',2)
                axis([0 4.5 0 4.5])

                subplot(2,2,2)
                title([num2str(ch) ': ' fcshdr.par(ch).name])

                XI_S = linspace(bhigh,prctile(datS_b,90),5);

                subplot(2,2,3)
                hold on
                plot(datS_b,datS_c,'.')
                YI_S = lsq_lut_piecewise(datS_b,datS_c,XI_S);
                plot(XI_S,YI_S,'ko-','LineWidth',2)
                axis([0 4.5 0 4.5])

                subplot(2,2,4)
                hold on
                plot(datS_b,datS_c_new,'r.')
                YI_S_new = lsq_lut_piecewise(datS_b,datS_c_new,XI_S);
                plot(XI_S,YI_S_new,'ko-','LineWidth',2)
                axis([0 4.5 0 4.5])
       
            end
        end
    end
    
    data = datS_new;


    %% compensate2
    
%     subplot(3,1,2)
%     
%     datSC = applyJCGate(SC_file,NT_file);
% %     datSC = datSC-repmat(autofluorescence,size(datSC,1),1);
%     
%     datSC_b = log10(datSC(:,10));
%     datSC_c = log10(datSC(:,11));
%     plot(datSC_b,datSC_c,'.')
%     
%     ind = datSC_b>3.2;
%     datSC_b = datSC_b(ind);
%     datSC_c = datSC_c(ind);
% %     plot(datSC_b,datSC_c,'.')
%     
%     p = polyfit(datSC_b,datSC_c,1)
%     s = datSC_b\datSC_c
%     
%     x = [-3:0.1:5];
%     hold on
%     plot(x, polyval(p,x), 'r-')
%     plot(x,x.*s,'r--')
%     plot((0-p(2))/p(1),0,'ro')
    
end

end