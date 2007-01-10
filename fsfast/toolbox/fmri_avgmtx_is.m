function Misa = fmri_avgmtx_is(Nch, Ns, XtXs)
%
% Misa = fmri_avgmtx_is(Nch, Ns)
% Misa = fmri_avgmtx_is(Nch, Ns, XtXs)
%
% Compute inter-session averaging matrix.
% Nch - number conditions * number of estimates per event
% Ns  - number of sessions
% XtXs  - list of (sum(x'x)) (one for each session).
%
% Misa will have dimension Nch x Nch x Ns
%
%


%
% fmri_avgmtx_is.m
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

if(nargin ~= 2 & nargin ~= 3)
  msg = 'Misa = fmri_avgmtx_is(Nch, Ns, <XtXs>)';
  qoe(msg);error(msg);
end

Misa = repmat(eye(Nch), [1 1 Ns]);


if(nargin == 2)
  Misa = Misa/Ns;
  return;
end

% Only reaches this point if nargin == 3 %

if(size(XtXs,1) ~= Nch | size(XtXs,3) ~= Ns)
  fprintf('size(XtXs) = '); fprintf(1,'%d  ',size(XtXs));fprintf(1,'\n');
  msg = sprintf('XtXs size is inconsistent with Nch=%d or Ns=%d',Nch,Ns);
  qoe(msg);error(msg);
end

dXtXssum = 0;
for s = 1:Ns,
  dXtXsum = dXtXsum + diag(XtXs(:,:,s));
end

for s = 1:Ns,
  d = diag(XtXs(:,:,s))./dXtXsum;
  Misa(:,:,s) = diag(d);
end

return;
