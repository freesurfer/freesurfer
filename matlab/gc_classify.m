function [datacid map p] = gc_classify(data,means,covar)
% [datacid map p] = gc_classify(data,means,covar)
% Apply gaussian classifier to a data set.
% datacid - estimates of the class identification
% map - maximum a posteriori p-value
% p - ndata-by-nclasses - prob that data is in given class
%  -- note: this is not a p-value. See gc_pvalue
% See also gc_synth_params, gc_synth_data, gc_train
%
% data  - nData x nVariates
% means - nClasses x nVariates
% covar - nVariates x nVariates x nClasses
%

[nClasses nVariates] = size(means);
nData = size(data,1);

% Scaling factor (not important except to interpret p)
e0 = (2*pi)^(nVariates/2);

p = zeros(nData,nClasses);
for nthClass = 1:nClasses
  dm = data - repmat(means(nthClass,:),[nData 1]);
  C = covar(:,:,nthClass);
  q = sum(dm .* (inv(C)*dm')',2); % Square of Mahalanobis Dist
  p(:,nthClass) = exp(-q/2)/(e0*sqrt(det(C)));
end

[map datacid] = max(p,[],2);

return;
