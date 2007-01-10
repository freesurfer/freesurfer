function [eres, sigest] = fmri_residual(y,X,hest)
%
% [eres sigest] = fmri_residual(y,X,hest)
%
% Computes the residual error between the actual signal y
% and the signal estimate sigest = X*h.
%
%


%
% fmri_residual.m
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

sigest = fmri_estsignal(X,hest);
eres = y - sigest;

return;
