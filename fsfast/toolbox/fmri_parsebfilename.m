function [stem, slice, ext] = fmri_getstem(bfilename)
% [stem slice ext] = fmri_getstem(bfilename)
% $Id: fmri_parsebfilename.m,v 1.1 2003/03/04 20:47:40 greve Exp $

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
