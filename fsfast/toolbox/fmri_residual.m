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

sigest = fmri_estsignal(X,hest);
eres = y - sigest;

return;
