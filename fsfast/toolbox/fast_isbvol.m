function isbvol = fast_isbvol(volid)
% isbvol = fast_isbvol(volid)
%
% Returns 1 if volume is in bfile format, ie, there exist file
% with name volid_%03d.bshort or .bfloat.  Actually, it just
% looks for the header volid_%03d.hdr.
% 
%


%
% fast_isbvol.m
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
