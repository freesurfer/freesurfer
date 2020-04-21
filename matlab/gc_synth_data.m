function [data datacid] = gc_synth_data(means,covar,nData)
% [data datacid] = gc_synth_data(means,covar,nData)
% Synthesize data to test a gaussian classifier. 
%
%  nData can be one of three things:
%     1. A single number indicating the nData per class
%     2. nClass numbers indicating the nData per class
%     3. ClassId list (nDataTotal = length(nData))
%  [means covar] = gc_synth_params(nClasses,nVariates)
%
% datacid = class membership
%
% See also gc_synth_params
% 

data = [];
if(nargin ~= 3)
  fprintf('data = gc_synth_data(means,covar,nData)\n');
  return;
end

[nClasses nVariates] = size(means);

resort = 0;
len_nData = length(nData);
if(len_nData == 1)
  % A single value
  nData = repmat(nData(:),[nClasses 1]);
elseif(len_nData ~= nClasses)
  % List of class for each data point
  cids = nData;
  cidlist = unique(cids); % Must be nClasses long
  for k=1:nClasses
    nData(k) = length(find(cids==cidlist(k)));
  end
  resort = 1;
end

data = [];
datacid = [];
for nthClass = 1:nClasses
  m = means(nthClass,:);
  C = covar(:,:,nthClass);
  F = chol(C)';
  r = F*randn(nVariates,nData(nthClass));
  r = r + repmat(m',[1 nData(nthClass)]);
  data = [data; r'];
  datacid = [datacid; nthClass*ones(nData(nthClass),1)];
end

% Resort if nData is class list
if(resort)
  data2    = zeros(size(data));
  datacid2 = zeros(size(datacid));
  for nthClass = 1:nClasses
    indA = find(datacid==cidlist(nthClass));
    indB = find(cids==cidlist(nthClass));
    data2(indB,:) = data(indA,:);
    datacid2(indB) = datacid(indA);
  end
  data = data2;
  datacid = datacid2;
end

return;

