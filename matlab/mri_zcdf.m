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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
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


p = (1 + erf(Z/sqrt(2)) )/2;

return;




