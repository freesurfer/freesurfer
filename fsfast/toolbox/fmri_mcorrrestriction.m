function RMC = fmri_mcorrestriction(RM,TR,nHEst,Delta,Tau)
%
% RMC = fmri_mcorrestriction(RM,TR,nHEst,Delta,Tau)
%
% Creates a restriction matrix for correlating with a
% hemodynmic response function (HRF).  The HRF is
% computed using a Gamma function with parameters 
% Delta and Tau over the estimation time window.
%
% RM is a restriction matrix created by CreateRM.
%
% See also: CreateRM(), HemoDyn()
%
%


%
% fmri_mcorrrestriction.m
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


t = TR*[0:nHEst-1]';
hHDIR = fmri_hemodyn(t,Delta,Tau);
nCond = size(RM,2)/nHEst; % excluding fix %
q = repmat(hHDIR',size(RM,1),nCond);
RMC = RM .* q;

return;
