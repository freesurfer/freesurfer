function fwhm = fast_ar12fwhm(ar1,d)
% fwhm = fast_ar12fwhm(ar1,d)
%
% Converts an AR1 to FWHM. 
% d is the voxel size. FWHM will be in units of d.
%
% $Id: fast_ar12fwhm.m,v 1.2 2007/02/02 05:10:56 greve Exp $

% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/02/02 05:10:56 $
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

fwhm = [];
if(nargin ~= 1 & nargin ~= 2)
  fprintf('ar1 = fast_ar12fwhm(ar1,d)\n');
  return;
end

fwhm = sqrt(log(256.0))*d./sqrt(-4*log(ar1));

return
