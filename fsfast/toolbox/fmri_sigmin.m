function [minsig, iminsig, indminsig] = fmri_minsig(sig)
% [minsig, iminsig, indminsig] = fmri_minsig(sig)
% 
% sig: - raw significance values (nrows X ncols X nplanes)


%
% fmri_sigmin.m
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

if(nargin ~= 1)
  msg = '[minsig, iminsig, indminsig] = fmri_minsig(sig)'
  qoe(msg); error(msg);
end

[minsig iminsig] = min(abs(sig),[],3);

%% Go through some gymnastics to get the signed pmin %%
[nrows ncols nplanes] = size(psig);
xx = [1:nrows]' * ones(1,ncols); %'
I1 = reshape1d(xx);
I2 = reshape1d(xx'); %'
indminsig = sub2ind(size(pSig),I1,I2,reshape1d(iminsig));
minsig    = reshape(sig(indminsig), [nrows ncols]);      

% Bonferoni Correction
minsig = minsig * size(sig,3);  

return;
