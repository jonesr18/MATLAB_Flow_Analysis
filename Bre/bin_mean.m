function [XI0, YI] = bin_mean(x,y,XI)


% vector for y measurements
YI = []; 
XI0 = [];

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
    yi = mean(y_);
    % take average of x edges
    xi = mean([XI(j-1) XI(j)]);
    
    % concatenate y measurements of bin
    YI = [YI; yi];
    XI0 = [XI0; xi];
end

YI = real(YI);

end