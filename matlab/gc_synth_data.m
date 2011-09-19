function [data datacid] = gc_synth_data(means,covar,nData)
% [data datacid] = gc_synth_data(means,covar,nData)
%   nData = nData per class.
%   [means covar] = gc_synth_params(nClasses,nVariates)
% datacid = class membership
% Synthesize data to test a gaussian classifier. 
% See also gc_synth_params


data = [];
if(nargin ~= 3)
  fprintf('data = gc_synth_data(means,covar,nData)\n');
  return;
end

[nClasses nVariates] = size(means);

data = [];
datacid = [];
for nthClass = 1:nClasses
  m = means(nthClass,:);
  C = covar(:,:,nthClass);
  F = chol(C)';
  r = F*randn(nVariates,nData);
  %tmp = r*r'/nData;
  r = r + repmat(m',[1 nData]);
  data = [data; r'];
  datacid = [datacid; nthClass*ones(nData,1)];
end

return;

