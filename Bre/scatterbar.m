function scatterbar(xdata, ydata, col, color)

figure(gcf)
hold on

dat = ydata;

datX = xdata;

mF = mean(datX);
r = mean([max(datX)-mF, mF-min(datX)]);
mr = 0.5/r;
        
densityplot(col+(datX-mF).*mr,lin2logicle(dat),color); 



% function scatterbar(Sfile, NTfile, ch, col, color)
% 
% figure(gcf)
% hold on
% 
% [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs(NTfile);
% dat = applyJCGate(Sfile,NTfile);
% 
% chF = getChannel(fcshdr,'FSC-A');
% datF = dat(:,chF);
% 
% mF = mean(datF);
% r = mean([max(datF)-mF, mF-min(datF)]);
% mr = 0.5/r;
%         
% densityplot(col+(datF-mF).*mr,lin2logicle(dat(:,ch)),color); 
% 
% 
% end
end