function [nelements,centers]=biexhist(Y,M)
%BIEXHIST Logiclly-scaled (Parks,et al.) histogram. 
%   BIEXHIST(Y) bins the elements of Y into 10 equally spaced containers
%   and produces a histogram plot of the results.  If Y is a
%   matrix, hist works down the columns.
% 
%   N = hist(Y,M), where M is a scalar, uses M bins.
%
%   Example:
%       [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs('sample.fcs');
%       GreenChannel = getChannel(fcshdr,'FIT');
%       GreenData = fcsdat(:,GreenChannel);
%       biexhist(greenData)
%
%   Written by
%   Breanna Stillo
%   bstillo@mit.edu
%   Last Updated: 2014-10-31;
%
%   UPDATES:
%   10/31 -- changed ylimits to reflect FlowJo.  Change default number of
%   bins


    ylog=lin2logicle(Y);
    
    if ~exist('M','var')
        spread=max(ylog)-min(ylog);
        M=round(77*spread);
    end
    [nelements,centers]=hist(ylog,M);
    maxN=round(max(nelements)*10^(-3))/(10^(-3));
    if maxN<max(nelements)
        maxN=maxN+1000;
    end
%     nelements=nelements./max(nelements).*100;
    h=area(centers,nelements,'FaceColor', [0.5 0.5 0.5],'LineWidth',1.5);
    
    TickVals=sort([-1e3 -[1:9].*10^2 -[1:9].*10^1 -[1:9].*10^0 0 [1:9].*10^0 [1:9].*10^1 [1:9].*10^2 [1:9].*10^3 [1:9].*10^4 1e5]);
    
    hha1 = gca;
    LargeTickVals=lin2logicle([-1e3 -1e2 -1e1 -1e0 0 1e0 1e1 1e2 1e3 1e4 1e5 ]);
    set(hha1, ...
       'ticklength',[0.02 0.09], ...
       'xlim',lin2logicle([-1500 262144]),...
       'ylim',[0 maxN],...%%%
       'xtick',LargeTickVals, ...
       'xticklabel',{},...
       'ytick',0:1000:maxN, ...
       'TickDir','out',...
       'LineWidth',1.5,'box','off')
   
   dist=-max(nelements)/15;
   text(lin2logicle(-1e3),dist,'-10^3','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14 )
   text(lin2logicle(0),dist,'0','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
   text(lin2logicle(1e3),dist,'10^3','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
   text(lin2logicle(1e4),dist,'10^4','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
   text(lin2logicle(1e5),dist,'10^5','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
   
    hha3=axes;
    set(hha3,'color','none', ...
        'XTick',lin2logicle(TickVals),'XTickLabel',{},...
        'YTick',0:1000:maxN,'YTickLabel',{},...
        'xlim',lin2logicle([-1500 262144]),'TickDir','out',...
    	'LineWidth',1)
    
    hha4=axes;
    set(hha4,'color','none', ...
        'XTick',[],'XTickLabel',{},...
        'YTick',[],'YTickLabel',{},...
        'xlim',lin2logicle([-1500 262144]),...
    	'box','on','LineWidth',1.5)
    
    linkaxes([hha1 hha3 hha4])
    
end