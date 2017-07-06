function biexaxis(axis,a)

%
%   Written by
%   Breanna Stillo
%   bstillo@mit.edu
%   Last Updated: 2015-10-17;


%         fprintf('Axis correction... may take a little time...please wait.\n')
%         hha1 = axis;
%         LargeTickVals=lin2logicle([-1e3 -1e2 -1e1 -1e0 0 1e0 1e1 1e2 1e3 1e4 1e5 ]);
%         set(hha1, ...
%            'ticklength',[0.02 0.09], ...
%            'xlim',lin2logicle([-1500 262144]),...
%            'xtick',LargeTickVals, ...
%            'xticklabel',{},...
%            'TickDir','out',...
%            'LineWidth',1.5,'box','off')
%        if ~exist('a','var')
%            set(hha1, ...
%            'ylim',lin2logicle([-1500 262144]),...
%            'ytick',LargeTickVals,...
%            'yticklabel',{})
%        end
%        
%        ylimvals = ylim;
%        if strcmp(get(axis,'yscale'), 'log')
%            dist=10^(log10(ylimvals(1))-(log10(ylimvals(2))-log10(ylimvals(1)))/18);
%        else
%            dist=ylimvals(1)-(ylimvals(2)-ylimvals(1))/18;
%        end
%        %x-axis
%        text(lin2logicle(-1e3),dist,'-10^3','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14 )
%        text(lin2logicle(0),dist,'0','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
%        text(lin2logicle(1e3),dist,'10^3','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
%        text(lin2logicle(1e4),dist,'10^4','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
%        text(lin2logicle(1e5),dist,'10^5','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
% 
%        xlimvals = xlim;
%        dist=xlimvals(1)-(xlimvals(2)-xlimvals(1))/18;
%        %y-axis
%        if ~exist('a','var')
%            text(dist,lin2logicle(-1e3),'-10^3','HorizontalAlignment','center','VerticalAlignment','middle' ,'FontSize',14)
%            text(dist,lin2logicle(0),'0','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
%            text(dist,lin2logicle(1e3),'10^3','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
%            text(dist,lin2logicle(1e4),'10^4','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
%            text(dist,lin2logicle(1e5),'10^5','HorizontalAlignment','center','VerticalAlignment','middle','FontSize',14)
%        end
%        
%         hha2=axes;
%         TickVals=sort([-1e3 -[1:9].*10^2 -[1:9].*10^1 -[1:9].*10^0 0 [1:9].*10^0 [1:9].*10^1 [1:9].*10^2 [1:9].*10^3 [1:9].*10^4 1e5]);
%         set(hha2,'color','none', ...
%             'XTick',lin2logicle(TickVals),'XTickLabel',{},...
%             'xlim',lin2logicle([-1500 262144]),...
%             'YTick',[],'YTickLabel',{},...
%             'TickDir','out',...
%             'LineWidth',1)
%         if ~exist('a','var')
%             set(hha2, ...
%             'YTick',lin2logicle(TickVals),'YTickLabel',{},...
%             'ylim',lin2logicle([-1500 262144]))
%         end
% 
%         hha3=axes;
%         set(hha3,'color','none', ...
%             'XTick',[],'XTickLabel',{},...
%             'YTick',[],'YTickLabel',{},...
%             'xlim',lin2logicle([-1500 262144]),...
%             'box','on','LineWidth',1.5)
%         if ~exist('a','var')
%             set(hha3,...
%             'ylim',lin2logicle([-1500 262144]))
%         end

    T=262144;   % The top data value
    M=4.5;      % Breadth of the display in decades
    r=-100;    % Negative range reference value

    W=(M-log10(T/abs(r)))/2;

    A=0;

%         fprintf('Axis correction... may take a little time...please wait.\n')
        fs = 10;
        hha1 = gca;
        p1 = get(gca, 'Position');
        LargeTickVals=lin2logicle([-1e3 -1e2 -1e1 -1e0 0 1e0 1e1 1e2 1e3 1e4 1e5 ]);
        set(hha1, ...
           'ticklength',[0.02 0.09], ...
           'xlim',lin2logicle([r T]),...%'ylim',lin2logicle([r T]),...
           'xtick',LargeTickVals, ...%'ytick',LargeTickVals,...
           'xticklabel',{'','','','','0','','','','10^{3}','10^{4}','10^{5}'},...%'yticklabel',{},
           'TickDir','out',...
           'LineWidth',1.5,'box','off')

        hha2=axes;
        TickVals=sort([-1e3 -[1:9].*10^2 -[1:9].*10^1 -[1:9].*10^0 0 [1:9].*10^0 [1:9].*10^1 [1:9].*10^2 [1:9].*10^3 [1:9].*10^4 1e5]);
        set(hha2,'color','none', ...
            'XTick',lin2logicle(TickVals),'XTickLabel',{},...
            'YTick',[],'YTickLabel',{},...
            'xlim',lin2logicle([r T]),'LineWidth',0.75,...
            'TickDir','out',...%'ylim',lin2logicle([r T]),'LineWidth',1,...
            'Position',p1)
        
        hha4 = axes('Position',get(hha2,'Position'),'box','on','xtick',[],'ytick',[],'LineWidth',1.5,'color','none');

        hha3=axes;
        set(hha3,'color','none', ...
            'XTick',[],'XTickLabel',{},...
            'YTick',[],'YTickLabel',{},...
            'xlim',lin2logicle([r T]),... %'ylim',lin2logicle([r T]),...
            'box','on','LineWidth',1.5,...
            'Position',p1)
        
        axes(hha1);
        
        linkaxes([hha1 hha2 hha3 hha4])
        
end