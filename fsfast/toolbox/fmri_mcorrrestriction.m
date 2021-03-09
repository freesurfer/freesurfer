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
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


t = TR*[0:nHEst-1]';
hHDIR = fmri_hemodyn(t,Delta,Tau);
nCond = size(RM,2)/nHEst; % excluding fix %
q = repmat(hHDIR',size(RM,1),nCond);
RMC = RM .* q;

return;
