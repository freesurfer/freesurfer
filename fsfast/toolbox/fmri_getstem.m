function [stem, ext] = fmri_getstem(bfilename)
% [stem ext] = fmri_getstem(bfilename)

if(nargin ~= 1)
  msg = 'USAGE: [stem ext] = fmri_getstem(bfilename)';
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
  stem = BFileName(1:ks-5);
else               
  ext = 'bfloat';
  stem = BFileName(1:kf-5);
end

return;
