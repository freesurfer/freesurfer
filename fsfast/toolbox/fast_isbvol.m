function isbvol = fast_isbvol(volid)
% isbvol = fast_isbvol(volid)
%
% Returns 1 if volume is in bfile format, ie, there exist file
% with name volid_%03d.bshort or .bfloat.  Actually, it just
% looks for the header volid_%03d.hdr.
% 
% $Id: fast_isbvol.m,v 1.1 2003/03/04 20:47:38 greve Exp $

isbvol = 0;

if(nargin ~= 1)
  msg = 'USAGE: r = fast_isbvol(volid)'
  qoe(msg); error(msg);
end

stem = deblank(volid);

nslices = 0;
for slice = 0:30
  fname = sprintf('%s_%03d.hdr',stem,nslices);
  fid = fopen(fname,'r');
  if(fid ~= -1) 
    fclose(fid);
    isbvol = 1;
    return;
  end
end

return;
