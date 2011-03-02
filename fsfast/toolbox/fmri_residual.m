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
%    $Date: 2011/03/02 00:04:06 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

sigest = fmri_estsignal(X,hest);
eres = y - sigest;

return;
