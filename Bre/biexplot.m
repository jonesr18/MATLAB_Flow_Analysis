function biexplot(x,y,S,F)

addpath(genpath('/Users/brestillo/Dropbox (MIT)/Weiss Lab/Work/Flow Analysis Repository MATLAB Package/release2 working'))

%BIEXPLOT(...) is the same as PLOT(...) except that logicle scales (Parks,
%et al.) are used for both the X- and Y- axes.
%   As with PLOT(...), various line types, plot symbols and colors may be obtained with
%   BIEXPLOT(X,Y,S,F) where S is a character string made from one element
%   from any or all the following 3 columns:
%
%          b     blue          .     point              -     solid
%          g     green         o     circle             :     dotted
%          r     red           x     x-mark             -.    dashdot 
%          c     cyan          +     plus               --    dashed   
%          m     magenta       *     star             (none)  no line
%          y     yellow        s     square
%          k     black         d     diamond
%          w     white         v     triangle (down)
%                              ^     triangle (up)
%                              <     triangle (left)
%                              >     triangle (right)
%                              p     pentagram
%                              h     hexagram
%
%   When the string 'density' is supplied for S, a density plot is created
%   using logicle scales.  This most closely imitates FlowJo plots
%
%   When F is set to 1, biexponential axis ticks are created
%
%   Example:
%       [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs('sample.fcs');
%       GreenChannel = getChannel(fcshdr,'FIT');
%       RedChannel = getChannel(fcshdr,'Red');
%       GreenData = fcsdat(:,GreenChannel);
%       RedData = fcsdat(:,RedChannel);
%       biexplot(greenData,redData,'density'
%
%   Written by
%   Breanna Stillo
%   bstillo@mit.edu
%   Last Updated: 2014-10-14;


    T=262144;   % The top data value
    M=4.5;      % Breadth of the display in decades
    r=-100;    % Negative range reference value

    W=(M-log10(T/abs(r)))/2;

    A=0;

    xlog=M.*logicleTransform(x,T,W,M,A);
    ylog=M.*logicleTransform(y,T,W,M,A);
%     xlog=lin2logicle(x);
%     ylog=lin2logicle(y);
    
    if ~exist('S','var')
        S='.';
        plot(xlog,ylog,S,'MarkerSize',2)
    else
        if strcmp(S,'density')
            densityplot(xlog,ylog)
        else
            plot(xlog,ylog,S,'MarkerSize',2)
        end
    end
    
    if exist('F','var')
    if F==1
        fprintf('Axis correction... may take a little time...please wait.\n')
        fs = 9;
        hha1 = gca;
        p1 = get(gca, 'Position');
        LargeTickVals=lin2logicle([-1e3 -1e2 -1e1 -1e0 0 1e0 1e1 1e2 1e3 1e4 1e5 ]);
        set(hha1, ...
           'ticklength',[0.02 0.09], ...
           'xlim',lin2logicle([r T]),...
           'ylim',lin2logicle([r T]),...
           'xtick',LargeTickVals, ...
           'ytick',LargeTickVals,...
           'xticklabel',{},...
           'yticklabel',{},'TickDir','out',...
           'LineWidth',1.5,'box','off')
       ylimvals = ylim;
       dist=ylimvals(1)-(ylimvals(2)-ylimvals(1))/36;
%        text(lin2logicle(-1e3),dist,'-10^3','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs )
       text(lin2logicle(0),dist,'0','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs)
       text(lin2logicle(1e3),dist,'10^3','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs)
       text(lin2logicle(1e4),dist,'10^4','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs)
       text(lin2logicle(1e5),dist,'10^5','HorizontalAlignment','center','VerticalAlignment','top','FontSize',fs)

%        text(dist,lin2logicle(-1e3),'-10^3','HorizontalAlignment','right','VerticalAlignment','middle' ,'FontSize',fs)
       text(dist,lin2logicle(0),'0','HorizontalAlignment','right','VerticalAlignment','middle','FontSize',fs)
       text(dist,lin2logicle(1e3),'10^3','HorizontalAlignment','right','VerticalAlignment','middle','FontSize',fs)
       text(dist,lin2logicle(1e4),'10^4','HorizontalAlignment','right','VerticalAlignment','middle','FontSize',fs)
       text(dist,lin2logicle(1e5),'10^5','HorizontalAlignment','right','VerticalAlignment','middle','FontSize',fs)

        hha2=axes;
        TickVals=sort([-1e3 -[1:9].*10^2 -[1:9].*10^1 -[1:9].*10^0 0 [1:9].*10^0 [1:9].*10^1 [1:9].*10^2 [1:9].*10^3 [1:9].*10^4 1e5]);
        set(hha2,'color','none', ...
            'XTick',lin2logicle(TickVals),'XTickLabel',{},...
            'YTick',lin2logicle(TickVals),'YTickLabel',{},...
            'xlim',lin2logicle([r T]),'TickDir','out',...
            'ylim',lin2logicle([r T]),'LineWidth',1,...
            'Position',p1)

        hha3=axes;
        set(hha3,'color','none', ...
            'XTick',[],'XTickLabel',{},...
            'YTick',[],'YTickLabel',{},...
            'xlim',lin2logicle([r T]),...
            'ylim',lin2logicle([r T]),'box','on','LineWidth',1.5,...
            'Position',p1)
        
        axes(hha1);
            
    end
    end
end