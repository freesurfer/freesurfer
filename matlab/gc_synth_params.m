function [means covar] = gc_synth_params(nClasses,nVariates)
% [means covar] = gc_synth_params(nClasses,nVariates)
% Synthesize parameters that can be used to synthesize data to test
% a gaussian classifier. 
% See also gc_synth_data, gc_train, gc_classify
%

if(nargin ~= 2)
  fprintf('[means covar] = gc_synth_params(nClasses,nVariates)\n');
  return
end

means = rand(nClasses,nVariates);
covar = zeros(nVariates,nVariates,nClasses);

for nthClass = 1:nClasses
  A = randn(nVariates,nVariates);
  covar(:,:,nthClass) = A'*A; % must be pos def
end

return;




