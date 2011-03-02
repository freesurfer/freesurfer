function err = fast_svcor(cor, cordir)
%
% err = fast_svcor(cor, cordir)
%


%
% fast_svcor.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:05 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
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
