function [p kappa F] = gc_pvalue(data,means,covar)
% [p kappa F] = gc_pvalue(data,means,covar)
% Compute p-value of data based on gaussian classifier. Note that
% this is different than the p computed by gc_classify.m. This
% p-value is the probably of seeing the data under the null
% of a class mean and covariance. kappa is the mahalanobis dist
%
% data  - nData x nVariates
% means - nClasses x nVariates
% covar - nVariates x nVariates x nClasses
%
% See also gc_synth_params, gc_synth_data, gc_train
%

[nClasses nVariates] = size(means);
nData = size(data,1);

p = zeros(nData,nClasses);
for nthClass = 1:nClasses
  dm = data - repmat(means(nthClass,:),[nData 1]);
  C = covar(:,:,nthClass);
  kappa(:,nthClass) = sqrt(sum(dm .* (inv(C)*dm')',2)); % mahalanobis distance
  F(:,nthClass) = (kappa(:,nthClass).^2)/nVariates;
  p(:,nthClass) = FTest(nVariates,1000,F(:,nthClass));
end

return;

