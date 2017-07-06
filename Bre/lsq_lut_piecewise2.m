function [ YI ] = lsq_lut_piecewise2( x, y, XI )
% LSQ_LUT_PIECEWISE Piecewise linear interpolation for 1-D interpolation (table lookup)
%   YI = lsq_lut_piecewise( x, y, XI ) obtain optimal (least-square sense)
%   vector to be used with linear interpolation routine.
%   The target is finding Y given X the minimization of function 
%           f = |y-interp1(XI,YI,x)|^2
%   
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
% A = sparse([]); 
A=[];

% vector for y measurements
Y = []; 

for j=2:length(XI)
    
    % get index of points in bin [XI(j-1) XI(j)]
    ix = x>=XI(j-1) & x<XI(j) & ~isnan(x) & ~isnan(y);
    
    % check if we have data points in bin
    if ~any(ix)
        warning(sprintf('Bin [%f %f] has no data points, check estimation. Please re-define X vector accordingly.',XI(j-1),XI(j)));
    end
    
    % get x and y data subset
    x_ = x(ix);
    y_ = y(ix);
    
    X = [ones(length(x_),1) x_];
    
    xnan = sum(isnan(x_))
    ynan = sum(isnan(y_))
    
    b=X\y_
    p = polyfit(x_,y_,1)
    
end
YI=0;
    
%     % create temporary matrix to be added to A
%     tmp = [(( -x_+XI(j-1) ) / ( XI(j)-XI(j-1) ) + 1) (( x_-XI(j-1) ) / ( XI(j)-XI(j-1) ))];
%     
%     % build matrix of measurement with constraints
%     [m1,n1]=size(A);
%     [m2,n2]=size(tmp);
%     A = [[A zeros(m1,n2-1)];[zeros(m2,n1-1) tmp]];
%     
%     % concatenate y measurements of bin
%     Y = [Y; y_];
% end
% 
% % obtain least-squares Y estimation
% % YI=A\Y;
% YI=mldivide(full(A),full(Y))


