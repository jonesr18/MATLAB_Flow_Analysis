function densityplot(x,y,base_col)
%DENSITYPLOT(X,Y) plots the vector Y vs vector X in dot-plot form where
%warmer colors indicate point density.
%   DENSITYPLOT(X,Y) creates a density plot of Y vs X where
%    X - is a vector of values
%    Y - is a vector of values the same size as X
%
%   Example:
%       [fcsdat, fcshdr, fcsdatscaled, fcsdatcomp] = fca_readfcs('sample.fcs');
%       GreenChannel = getChannel(fcshdr,'FIT');
%       RedChannel = getChannel(fcshdr,'Red');
%       GreenData = fcsdat(:,GreenChannel);
%       RedData = fcsdat(:,RedChannel);
%       densityplot(greenData,redData);
%
%   Written by
%   Breanna Stillo
%   bstillo@mit.edu
%   Last Updated: 2014-10-14
    

    fprintf('Density plot generation... may take a little time... please wait.\n')
    xdist=max(x)-min(x);
    ydist=max(y)-min(y);

    xdens=length(x)/xdist;
    ydens=length(x)/ydist;
    dX=100/xdens;
    dY=100/ydens;
    
    
    neighbors=zeros(length(x),1);
    for j=1:length(x)
        xval=x(j);
        yval=y(j);

        neighbors(j)=length(find(x>(xval-dX) & x<(xval+dX) & y>(yval-dY) & y<(yval+dY)));

    end
    
    highN=prctile(neighbors,90);
    
    minN=min(neighbors);
    
    col=zeros(length(x),3);
    for j=1:length(x)
        if ~exist('base_col','var')
            col(j,:)=getColor(neighbors(j),minN,highN);
        else
            col(j,:)=getColor(neighbors(j),minN,highN,base_col);
        end
    end
    scatter(x,y,3,col,'filled')


end


function c= getColor(v,vmin,vmax,base_col)

       if (v < vmin)
          v = vmin;
       elseif (v > vmax)
          v = vmax;
       end
       dv = vmax - vmin;

       
    if ~exist('base_col','var')
       c = [1,1,1];

       if (v < (vmin + 0.25 * dv)) 
          c(1) = 0;
          c(2) = 4 * (v - vmin) / dv;
       elseif (v < (vmin + 0.5 * dv)) 
          c(1) = 0;
          c(3) = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
       elseif (v < (vmin + 0.75 * dv)) 
          c(1) = 4 * (v - vmin - 0.5 * dv) / dv;
          c(3) = 0;
       else
          c(2) = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
          c(3) = 0;
       end
    else
        c = base_col;
        dv = dv*1.5;
    
       if (v < (vmin + 0.25 * dv))  % only a few neighbors
          c = c - c .* (v - vmin) / dv;
       elseif (v < (vmin + 0.5 * dv)) % bunch of neighbors
          c = c - c .* (v - vmin) / dv;
       elseif (v < (vmin + 0.75 * dv)) % many neighbors
          c = c - c .* (v - vmin) / dv;
       else %so many neighbors
          c = c - c .* (v - vmin) / dv;
       end
    end

end