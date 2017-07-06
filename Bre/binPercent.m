function [ XInew, YI ] = binPercent( x, y, XI, P )

%   INPUT
%       x measured data vector
%       y measured data vector
%       XI break points of 1-D table
%
%   OUTPUT
%       YI interpolation points of 1-D table
%           y = interp1(XI,YI,x)
%


if size(x,2) ~= 1
    error('Vector x must have dimension n x 1.');   
elseif size(y,2) ~= 1
    error('Vector y must have dimension n x 1.');    
elseif size(x,1) ~= size(x,1)
    error('Vector x and y must have dimension n x 1.'); 
end

% matrix defined by x measurements
A = sparse([]); 

% vector for y measurements
XInew = []; 
YI = []; 

for j=2:length(XI)
    
    % get index of points in bin [XI(j-1) XI(j)]
    ix = x>=XI(j-1) & x<XI(j);
    
    % check if we have data points in bin
    if ~any(ix)
        warning(sprintf('Bin [%f %f] has no data points, check estimation. Please re-define X vector accordingly.',XI(j-1),XI(j)));
    end
    
    % get x and y data subset
    x_ = x(ix);
    y_ = y(ix);
       
    
    % concatenate y measurements of bin
    XInew = [XInew mean(x_)];
    YI = [YI sum(y_>P)/length(y_)];
end


