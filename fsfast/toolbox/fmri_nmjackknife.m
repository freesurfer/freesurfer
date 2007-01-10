function a = fmri_nmjackknife(v)
% Jackknife averaging with normalization
%
% a = fmri_nmjackknife(v)
%
% size(v) = (Nh,Nvox,Nsamples)


%
% fmri_nmjackknife.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
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

[Nh Nvox Nsamples] = size(v);

for s = 1:Nsamples
  jk = find([1:Nsamples] ~= s);
  vtemplate  = fmri_norm(mean(v(:,:,jk),3),2);
  %vtemplate  = fmri_norm(randn(size(v(:,:,1))),2);

  if(size(v,1) > 1)
    a(s,:) = sum(v(:,:,s).*vtemplate);
  else
    a(s,:) = v(:,:,s).*vtemplate;
  end

end

return;
