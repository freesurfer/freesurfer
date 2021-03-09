function p = mri_zcdf(Z)
% p = mri_zcdf(Z)
% Cumulative distribution function for z, ie,
%   p = prob that a number drawn from a z distribution is <= Z.
%
% Simply returns:
%   p = (1 + erf(Z/sqrt(2)) )/2;
%
% To use this as a z-to-p converter, compute
%   p = 1 - mri_zcdf(z);
% 


%
% mri_zcdf.m
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


p = (1 + erf(Z/sqrt(2)) )/2;

return;




