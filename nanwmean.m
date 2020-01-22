function y = nanwmean(x,w,dim)
%NANWMEAN   Weighted Average or mean value.
%   For vectors, NANWMEAN(X,W) is the weighted mean value of the elements in X
%   using non-negative weights W. For matrices, NANWMEAN(X,W) is a row vector 
%   containing the weighted mean value of each column.  For N-D arrays, 
%   NANWMEAN(X,W) is the weighted mean value of the elements along the first 
%   non-singleton dimension of X.
%
%   Each element of X requires a corresponding weight, and hence the size 
%   of W must match that of X.
%
%   NANWMEAN(X,W,DIM) takes the weighted mean along the dimension DIM of X. 
%
%   Class support for inputs X and W:
%      float: double, single
%
%   Example:
%       x = rand(5,2);
%       w = rand(5,2);
%       wmean(x,w)

if nargin<2
    error('Not enough input arguments.');
end

% Remove not finite entries:
ind_bad = find(~isfinite(x) | ~isfinite(w));
x(ind_bad) = [];
w(ind_bad) = [];

% If there are no good values:
if isempty(x) || isempty(w)
   y = 0;
   return
end

% Check that dimensions of X match those of W.
if(~isequal(size(x), size(w)))
    error('Inputs x and w must be the same size.');
end

% Check that all of W are non-negative.
if (any(w(:)<0))
    error('All weights, W, must be non-negative.');
end

% Check that there is at least one non-zero weight.
if (all(w(:)==0))
    error('At least one weight must be non-zero.');
end

if nargin==2, 
  % Determine which dimension SUM will use
  dim = min(find(size(x)~=1));
  if isempty(dim), dim = 1; end
end

y = sum(w.*x,dim)./sum(w,dim);