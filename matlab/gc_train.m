function [means covar] = gc_train(data,datacid);
% [means covar] = gc_synth_params(data,datacid);
% Estimates gaussian classifier parameters from a training set.
% See also gc_synth_params, gc_synth_data, gc_classify
%
% data  - nData x nVariates
% means - nClasses x nVariates
% covar - nVariates x nVariates x nClasses
%

cids = unique(datacid);
nClasses = length(cids);
nVariates = size(data,2);

clear means covar;
for nthClass = 1:nClasses
  ind = find(datacid == cids(nthClass));
  nDC = length(ind);
  y = data(ind,:);
  m = mean(y);
  means(nthClass,:) = m;
  r = y - repmat(m,[nDC 1]);
  C = (r'*r)/nDC;
  covar(:,:,nthClass) = C;
end

return;

