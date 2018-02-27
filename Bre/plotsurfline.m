function plotsurfline(varargin)
%PLOTSURFLINE  plots a 2d array as a srface plot with lines along either the x
%or y direction
%   PLOTSURFLINE(xvals, yvals, Z, linedir) plots the values in Z using the
%   two vector arguments xvals and yvals must have length(x) = n and 
%   length(y) = m where [n,m] = size(Z).
%
%   PLOTSURFLINE(Z)	Z is a 2D array of values arranged (x,y). uses xvals =
%   1:n and yvals =1:m.
%
%   PLOTSURFLINE(__,Name,Value)  specifies the plotsurfline properties using one or
%   more Name, Value pair arguments specified below.
%
%   PLOTSURFLINE(ax,__)  plots into the axes specified by ax instead of in
%   the current axes (gca). The option ax can precede any of the input
%   argument combinations in the previous syntaxes.
%
%
%   Properties: 'linedir'   :   {'x','y'}
%               'colormap'  :   {char} e.g. 'parula'

    
    switch nargin
        case 0
            % not allowed
            error('Not enough input arguments.')
        case 1
            % should just be Z
            if isa(varargin{1},'double')
                Z = varargin{1};
                ax = gca;
                xvals = 1:size(Z,1);
                yvals = 1:size(Z,2);
            else
                error('Not enough input arguments.')
            end
        case 2
            % should be ax followed by Z
            if isa(varargin{1},'matlab.graphics.axis.Axes') && isa(varargin{2},'double')
                ax = varargin{1};
                Z = varargin{2};
                xvals = 1:size(Z,1);
                yvals = 1:size(Z,2);
            else
                error('Error in argument')
            end
        case 3
            % either (xvals,yvals,Z) or (Z,Name,Value)
            if isvector(varargin{1}) && isvector(varargin{2}) && ismatrix(varargin{3}) &&...
                    isa(varargin{1},'double') && isa(varargin{2},'double') && isa(varargin{3},'double')
                xvals = varargin{1};
                yvals = varargin{2};
                Z = varargin{3};
                ax = gca;
            elseif ismatrix(varargin{1}) && isa(varargin{1},'double')
                Z = varargin{1};
                xvals = 1:size(Z,1);
                yvals = 1:size(Z,2);
                ax = gca;
                props = varargin(2:nargin);
            else
                error('Error in argument')
            end
        case 4
            %(ax,Z,Name,Value) or (ax,xvals,yvals,Z)
            if isa(varargin{1},'matlab.graphics.axis.Axes') && isvector(varargin{2})...
                    && isvector(varargin{3}) && ismatrix(varargin{4})
                ax = varargin{1};
                xvals = varargin{2};
                yvals = varargin{3};
                Z = varargin{4};
            elseif isa(varargin{1},'matlab.graphics.axis.Axes') && ismatrix(varargin{2})
                props = varargin(3:nargin);
                ax = varargin{1};
                Z = varargin{2};
                xvals = 1:size(Z,1);
                yvals = 1:size(Z,2);
            else
                error('Error in argument')
            end
        otherwise %nargs>4
            % (Z,Name,Value...) or (xvals, yvals, Z, Name, Value) or (ax,
            % Z, Name, Value...) or (ax, xvals, yvals, Z, Name, Value...)
            if isa(varargin{1},'matlab.graphics.axis.Axes')
                if isvector(varargin{2}) && isvector(varargin{3}) && ismatrix(varargin{4})
                    ax = varargin{1};
                    xvals = varargin{2};
                    yvals = varargin{3};
                    Z = varargin{4};
                    props = varargin(5:nargin);
                elseif ismatrix(varargin{2})
                    ax = varargin{1};
                    Z = varargin{2};
                    xvals = 1:size(Z,1);
                    yvals = 1:size(Z,2);
                    props = varargin(3:nargin);
                else
                    error('Error in argument')
                end
            elseif isvector(varargin{1}) && isvector(varargin{2}) && ismatrix(varargin{3})
                xvals = varargin{1};
                yvals = varargin{2};
                Z = varargin{3};
                ax = gca;
                props = varargin(4:nargin);
            elseif ismatrix(varargin{1}) && isa(varargin{1},'double')
                
                Z = varargin{1};
                xvals = 1:size(Z,1);
                yvals = 1:size(Z,2);
                ax = gca;
                props = varargin(2:nargin);
            else
                error('Error in argument')
            end
                
    end

    %check input dimensions
    if ~(length(xvals) == size(Z,1) && length(yvals) == size(Z,2))
        error('Dimension mismatch. xvals and yvals must have length(x) = n and length(y) = m where [n,m] = size(Z)')
    end
    
    %check if linedir  or cmap assigned
    if exist('props','var')
        for p = 1:2:numel(props)-1
            if ~isa(props{p},'char')
                error('Property must be a char')
            end
            switch props{p}
                case 'linedir'
                    linedir = props{p+1};
                case 'colormap'
                    cmap = props{p+1};
                otherwise
                    error('Property does not exist')
            end
        end
    end
    if ~exist('linedir','var')
        linedir = 'x';
    end
    if ~exist('cmap','var')
        cmap = 'cool';
    end

    
    if ismatrix(Z)
        hold on
        Nx = size(Z,1);
        Ny = size(Z,2);
        switch linedir
            case 'x'
                colorscale = eval([cmap '(' num2str(Ny) ')']);
                surf(xvals,yvals,Z','FaceColor',[0.9 0.9 0.9],'EdgeColor','none','facealpha',0.9)
                for xi = 1:Nx
                    plot3(ones(Ny,1).*xvals(xi),yvals,Z(xi,:),'k-','Linewidth',0.5)
                end
                for yi = 1:Ny
                    plot3(xvals,ones(Nx,1).*yvals(yi),Z(:,yi),'.-','color',colorscale(yi,:),'Markersize',15,'Linewidth',2)
                end
            case 'y'
                colorscale = eval([cmap '(' num2str(Nx) ')']);
                surf(xvals,yvals,Z','FaceColor',[0.9 0.9 0.9],'EdgeColor','none','facealpha',0.9)
                for yi = 1:Ny
                    plot3(xvals,ones(Nx,1).*yvals(yi),Z(:,yi),'k-','Linewidth',0.5)
                end
                for xi = 1:Nx
                    plot3(ones(Ny,1).*xvals(xi),yvals,Z(xi,:),'.-','color',colorscale(xi,:),'Markersize',15,'Linewidth',2)
                end
                
        end
        view(-45,30)
    else
        error('This function can only plot 2 dimensional data')
    end
        


end