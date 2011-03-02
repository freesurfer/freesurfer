function [stem, slice, ext] = fmri_getstem(bfilename)
% [stem slice ext] = fmri_getstem(bfilename)
%


%
% fmri_parsebfilename.m
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
