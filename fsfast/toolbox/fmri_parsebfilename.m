function [stem, slice, ext] = fmri_getstem(bfilename)
% [stem slice ext] = fmri_getstem(bfilename)
%


%
% fmri_parsebfilename.m
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

if(nargin ~= 1)
  msg = 'USAGE: [stem slice ext] = fmri_getstem(bfilename)';
  qoe(msg); error(msg);
end

BFileName = deblank(bfilename);
ks = findstr(BFileName,'.bshort');
kf = findstr(BFileName,'.bfloat');

if(isempty(ks) & isempty(kf))
  msg = 'BFileName must be either bshort or bfloat';
  qoe(msg); error(msg);
end

if( ~isempty(ks) ) 
  ext = 'bshort';
  slice = str2num(BFileName(ks-3:ks));
  stem = BFileName(1:ks-5);
else               
  ext = 'bfloat';
  slice = str2num(BFileName(kf-3:kf));
  stem = BFileName(1:kf-5);
end

return;
