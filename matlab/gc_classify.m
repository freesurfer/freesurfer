function [datacid map p] = gc_classify(data,means,covar)
% [datacid map p] = gc_classify(data,means,covar)
% Apply gaussian classifier to a data set.
% datacid - estimates of the class identification
% map - maximum a posteriori p-value
% p - ndata-by-nclasses - p-value for each class
% See also gc_synth_params, gc_synth_data, gc_train
%
% $Id: gc_classify.m,v 1.2 2011/09/26 18:07:15 greve Exp $

[nClasses nVariates] = size(means);
nData = size(data,1);

% Scaling factor (not important except to interpret p)
e0 = (2*pi)^(nVariates/2);

p = zeros(nData,nClasses);
for nthClass = 1:nClasses
  dm = data - repmat(means(nthClass,:),[nData 1]);
  C = covar(:,:,nthClass);
  q = sum(dm .* (inv(C)*dm')',2);
  p(:,nthClass) = exp(-q/2)/(e0*sqrt(det(C)));
end

[map datacid] = max(p,[],2);

return;

% Old and slow
% p = zeros(nData,nClasses);
% for nthData = 1:nData
%   d = data(nthData,:);
%   for nthClass = 1:nClasses
%     dm = d - means(nthClass,:);
%     C = covar(:,:,nthClass);
%     q = dm*inv(C)*dm';
%     p(nthData,nthClass) = exp(-q/2)/(e0*sqrt(det(C)));
%   end
% end













