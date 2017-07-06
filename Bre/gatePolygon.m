function [indices, PolyVert] = gatePolygon(xdat,ydat,axisScales,position)
% gates a given population.  Uses the polygon vertices specified by
% 'position' or if this is not given, the user is asked to specify the polygon
% through UI

plotsOn=0;

if nargin<4  % positions are not given --> use user input
%     figure
    if strcmp(axisScales,'loglog')
        loglog(xdat,ydat,'.','MarkerSize',2)
        set(gca, 'XScale', 'log','YScale', 'log')
    elseif strcmp(axisScales,'semilogy')
        semilogy(xdat,ydat,'.','MarkerSize',2)
        set(gca, 'YScale', 'log')%,'xlim', [0 250000])
    elseif strcmp(axisScales,'semilogx')
        semilogx(xdat,ydat,'.','MarkerSize',2)
        set(gca, 'XScale', 'log')
    else
        plot(xdat,ydat,'.','MarkerSize',2)
    end
    hold on
    h = impoly;
    position = wait(h);
    plotsOn=1;
end

numVertex = length(position(:,1));

[minX minXI]=min(position(:,1));
if minXI == numVertex
    fwdI = 1;
else
    fwdI = minXI + 1;
end
if minXI == 1
    prevI = numVertex;
else
    prevI = minXI - 1;
end

if position(fwdI,2)>position(prevI,2)
    dir=1;
else
    dir=-1;
end


% trace top curve from left to right and find points under it

currI = minXI;
nextI = currI+dir;
if nextI == 0
    nextI = numVertex;
elseif nextI>numVertex
    nextI = 1;
end

xpH=[];

while position(nextI,1)>position(currI,1)
    
    if strcmp(axisScales,'loglog') || strcmp(axisScales,'semilogy')
        ysc='log';%get(gca,'YScale');
    else
        ysc = 'linear';
    end
    if strcmp(ysc,'linear')
        p=polyfit(position([currI,nextI],1),position([currI,nextI],2),1);
        yeval=polyval(p,xdat);  
        xpH=[xpH;find(xdat>position(currI,1) & xdat<position(nextI,1) & yeval-ydat>=0)];
    else
        p=polyfit(position([currI,nextI],1),log10(position([currI,nextI],2)),1);
        yeval=polyval(p,xdat);  
        xpH=[xpH;find(xdat>position(currI,1) & xdat<position(nextI,1) & yeval-log10(ydat)>=0)];
    end
    currI = nextI;
    nextI = currI+dir;
    if nextI == 0
        nextI = numVertex;
    elseif nextI>numVertex
        nextI = 1;
    end
end

% trace bottom curve from left to right and find points above it

currI = minXI;
nextI = currI-dir;
if nextI == 0
    nextI = numVertex;
elseif nextI>numVertex
    nextI = 1;
end

xpL=[];

while position(nextI,1)>position(currI,1)
    
    if strcmp(axisScales,'loglog') || strcmp(axisScales,'semilogy')
        ysc='log';%get(gca,'YScale');
    else
        ysc = 'linear';
    end
    if strcmp(ysc,'linear')
        p=polyfit(position([currI,nextI],1),position([currI,nextI],2),1);   
        yeval=polyval(p,xdat);
        xpL=[xpL;find(xdat>position(currI,1) & xdat<position(nextI,1) & yeval-ydat<=0)];
    else
        p=polyfit(position([currI,nextI],1),log10(position([currI,nextI],2)),1);   
        yeval=polyval(p,xdat);
        xpL=[xpL;find(xdat>position(currI,1) & xdat<position(nextI,1) & yeval-log10(ydat)<=0)];
    end
        
    currI = nextI;
    nextI = currI-dir;
    if nextI == 0
        nextI = numVertex;
    elseif nextI>numVertex
        nextI = 1;
    end
end

xpI=intersect(xpL,xpH);

if plotsOn
    hold on
    semilogy(xdat(xpI),ydat(xpI),'r.','MarkerSize',2)
    hold off
end

indices=xpI;


PolyVert=position;

end