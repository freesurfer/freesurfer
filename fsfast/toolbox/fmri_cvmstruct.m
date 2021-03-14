function cvmstruct = fmri_cvmstruct;
% cvmstruct = fmri_cvmstruct
% Creates covariance matrix structure


%
% fmri_cvmstruct.m
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

cvmstruct.version = 2;  % version of this header
cvmstruct.n       = 0;  % number of data points
cvmstruct.d       = 0;  % distance between data points (eg, TR)
cvmstruct.sz      = 0;  % size of cvm (rows or columns);
cvmstruct.cvm     = 0;  % actual covariance matrix
cvmstruct.norm    = 0;  % 1 = has been normalized
cvmstruct.inv     = 0;  % 1 = has been inverted
cvmstruct.acoravg = []; % average autocor function (with norm)
cvmstruct.acorstd = []; % stddev autocor function  (with norm)

return;
