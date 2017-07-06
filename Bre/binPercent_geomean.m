function [XI0, YI] = binPercent_geomean(x,y,N)


% vector for y measurements
YI = []; 
XI0 = [];

% if length(XI) == 1
    XI = prctile(x,linspace(0,100,N));
% end

for j=2:length(XI)
    
    % get index of points in bin [XI(j-1) XI(j)]
%     ix = x>=XI(j-1) & x<XI(j);
    ix = x>=XI(j-1) & x<XI(j) & ~isnan(x) & ~isnan(y);
    
    % check if we have data points in bin
    if ~any(ix)
        warning(sprintf('Bin [%f %f] has no data points, check estimation. Please re-define X vector accordingly.',XI(j-1),XI(j)));
    end
    
    % get x and y data subset
    x_ = x(ix);
    y_ = y(ix);
    
    % take geomean of y subset
    yi = geomean_neg(y_);
    % take average of x edges
    xi = geomean_neg(x_);
    
    % concatenate y measurements of bin
    YI = [YI; yi];
    XI0 = [XI0; xi];
end

YI = real(YI);

end