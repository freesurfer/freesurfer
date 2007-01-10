function scm = fmri_seq2scm(seq, Nh)
%
% scm = fmri_seq2scm(seq, Nh)
%
% Converts a sequence of stimuli to a stimulus convolution
% matrix assuming that the time between presentations equals
% the TR.
%
%


%
% fmri_seq2scm.m
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

[nseq nruns] = size(seq);

v0 = zeros(nseq,1);
r0 = zeros(1,Nh);

nNNCond = max(seq) - min(seq);

%scm = zeros(nseq,nNNCond*Nh,nruns);
scm = [];

for r = 1:nruns,

  x = [];
  for cond = min(seq(:,r))+1 : max(seq(:,r))
    ind = find(seq(:,r)==cond);
    c0  = v0;
    c0(ind) = 1;
    r0(1) = c0(1);
    xc = toeplitz(c0,r0);
    x = [x xc];
  end

  scm(:,:,r) = x;
end

return;
