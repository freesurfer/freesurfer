function err = fast_svcor(cor, cordir)
%
% err = fast_svcor(cor, cordir)
%

err = 1;

if(nargin ~= 2)
  fprintf('USAGE: err = fast_svcor(cor, cordir)\n');
  return;
end

ncorslices = size(cor,2);

for sliceno = 1:ncorslices
  fprintf('%3d ',sliceno);
  corslicefile = sprintf('%s/COR-%03d',cordir,sliceno);
  err = fast_svcorslice(squeeze(cor(:,sliceno,:)),corslicefile);
  if(err) return; end
end

fprintf('\n');

err = 0;

return;
