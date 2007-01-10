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
%    $Date: 2007/01/10 22:02:32 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
