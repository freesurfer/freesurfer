function p = gc_pvalue(data,means,covar)
% pvalue = gc_pvalue(data,means,covar)
% Compute p-value of data based on gaussian classifier. Note that
% this is different than the p computed by gc_classify.m. This
% p-value is the probably of seeing the data under the null
% of a class mean and covariance.
%
% data  - nData x nVariates
% means - nClasses x nVariates
% covar - nVariates x nVariates x nClasses
%
% See also gc_synth_params, gc_synth_data, gc_train
%
% $Id: gc_pvalue.m,v 1.1 2012/04/20 16:04:37 greve Exp $

[nClasses nVariates] = size(means);
nData = size(data,1);

p = zeros(nData,nClasses);
for nthClass = 1:nClasses
  dm = data - repmat(means(nthClass,:),[nData 1]);
  C = covar(:,:,nthClass);
  F = sum(dm .* (inv(C)*dm')',2)/nVariates;
  p(:,nthClass) = FTest(nVariates,1000,F);
end

return;

