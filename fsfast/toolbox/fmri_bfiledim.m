function [nrows, ncols, ntp, fs, ns, endian, bext] = fmri_bfiledim(stem)
% [nrows ncols ntp fs ns endian bext] = fmri_bfiledim(stem)
%
%


%
% fmri_bfiledim.m
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

if(nargin ~= 1) 
  msg = '[nrows ncols ntp fs ns endian bext] = fmri_bfiledim(stem)';
  qoe(msg); error(msg);
end

stem = deblank(stem);

fs = 0;
fname = sprintf('%s_%03d.hdr',stem,fs);
fid = fopen(fname,'r');
while(fid == -1 & fs < 100)
  fs = fs + 1;
  fname = sprintf('%s_%03d.hdr',stem,fs);
  fid = fopen(fname,'r');
end

if(fs == 100)
  msg = sprintf('Could not find any slices with stem %s',stem);
  qoe(msg); error(msg);
end

hdr = fscanf(fid,'%d',[1,4]);
fclose(fid);
nrows  = hdr(1);
ncols  = hdr(2);
ntp    = hdr(3);
endian = hdr(4);

ns = 1;
slice = fs + ns;
fname = sprintf('%s_%03d.hdr',stem,slice);
fid = fopen(fname,'r');
while(fid ~= -1)
  fclose(fid);
  ns = ns + 1;
  slice = fs + ns;
  fname = sprintf('%s_%03d.hdr',stem,slice);
  fid = fopen(fname,'r');
end


% Get extension %
bext = 'bshort';
fname = sprintf('%s_%03d.%s',stem,fs,bext);
fid = fopen(fname,'r');
if(fid == -1)
  bext = 'bfloat';
  fname = sprintf('%s_%03d.%s',stem,fs,bext);
  fid = fopen(fname,'r');
  if(fid == -1)
     msg = sprintf('Could not find bshort or bfloat for %s',stem);
     qoe(msg);error(msg);
  end
end

fclose(fid);

return;
