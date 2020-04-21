function X = fast_outliermtx(outliers,ndata)
% X = fast_outliermtx(outliers,ndata)
% 
% Constructs a design matrix to remove outliers
% outliers - list of 1-based outlier numbers
% ndata - number of data points (eg, time points)
%

X = [];
if(nargin ~= 2)
  fprintf('X = fast_outliermtx(outliers,ndata)\n');
  return;
end

noutliers = length(outliers);
if(noutliers == 0) return; end

if(max(outliers) > ndata)
  fprintf('ERROR: max(outliers) > ndata\n');
  return;
end

c = 1:noutliers;

X = zeros(ndata,noutliers);
ind = sub2ind(size(X),outliers(:),c(:));
X(ind) = 1;

return


  







