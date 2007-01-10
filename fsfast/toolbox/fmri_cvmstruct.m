function cvmstruct = fmri_cvmstruct;
% cvmstruct = fmri_cvmstruct
% Creates covariance matrix structure


%
% fmri_cvmstruct.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
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
