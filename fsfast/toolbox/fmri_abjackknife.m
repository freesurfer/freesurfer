function [a, vtemplate] = fmri_abjackknife(v,alist,blist)
% [a vtemplate] = fmri_abjackknife(v,alist,blist)
%
% size(v) = (Nh,Nvox,Nsamples)


%
% fmri_abjackknife.m
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

if(nargin ~= 3)
  msg = 'USAGE: a = fmri_abjackknife(v,alist,blist)';
  qoe(msg);error(msg);
end

[Nh Nvox Nsamples] = size(v);

vtemplate  = fmri_norm(mean(v(:,:,blist),3),2);

n = 1;
for s = alist,

  if(size(v,1) > 1)
    a(n,:) = sum(v(:,:,s).*vtemplate);
  else
    a(n,:) = v(:,:,s).*vtemplate;
  end
  n = n + 1;
end

return;
