function [nslices, nrows, ncols, nt, endian, bext, hdrdat] = fmri_bvoldim(stem)
% [nslices nrows ncols nt endian bext hdrdat] = fmri_bvoldim(stem)
% 
%


%
% fmri_bvoldim.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
%    $Revision: 1.3 $
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

nslices = 0;
nrows   = 0;
ncols   = 0;
nt      = 0;
endian  = -1;
bext    = '';
hdrdat  = [];

if(nargin ~= 1)
  msg = 'USAGE: [nslices nrows ncols nt] = fmri_bvoldim(stem)';
  qoe(msg); error(msg);
end

stem = deblank(stem);

firstslice = -1;
nslices = 0;
for slice = 0:1000
  fname = sprintf('%s_%03d.hdr',stem,nslices);
  fid = fopen(fname,'r');
  if(fid == -1) break;  end
  hdr = fscanf(fid,'%d',[1,4]);
  if(isempty(hdr))
    fprintf('ERROR: reading %s\n',fname);
    nslices = 0;
    return;
  end
  nrows = hdr(1);
  ncols = hdr(2);
  nt = hdr(3);
  endian = hdr(4);
  fclose(fid);

  % Get the extension %
  if(firstslice == -1) 
    firstslice = nslices ;
    bext = 'bshort';
    fname = sprintf('%s_%03d.%s',stem,nslices,bext);
    fid = fopen(fname,'r');
    if(fid == -1) 
      bext = 'bfloat';
      fname = sprintf('%s_%03d.%s',stem,nslices,bext);
      fid = fopen(fname,'r');
      if(fid == -1) 
        nslices = 0;
        nrows   = 0;
        ncols   = 0;
        nt      = 0;
        bext = '';
        endian = -1;
        return;
      end
    end
    fclose(fid);
  end

  nslices = nslices + 1;
end

fname = sprintf('%s.dat',stem);
fid = fopen(fname,'r');
if(fid ~= -1)
  fclose(fid);
  hdrdat = fmri_lddat3(fname);
end

return;
