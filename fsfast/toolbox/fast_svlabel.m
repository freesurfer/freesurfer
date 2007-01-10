function err = fast_svlabel(vtxno, xyz, vstat, hdr, filename)
% err = fast_svlabel(vtxno, xyz, vstat, hdr, filename)
%


%
% fast_svlabel.m
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

err = 0;

if(nargin ~= 5)
  msg = 'USAGE: err = fast_svlabel(vtxno, xyz, vstat, hdr, filename)';
  err = 1;
  qoe(msg);error(msg);
end

nvtxs = size(xyz,1);
if(isempty(vtxno)) vtxno = zeros(nvtxs,1); end
if(isempty(vstat)) vstat = zeros(nvtxs,1); end

if(nvtxs  ~= length(vtxno))
  msg = sprintf('nvtxs = %d, != length(vtxno) = %d\n',nvtxs,length(vtxno));
  err = 1;
  qoe(msg);error(msg);
end

if(nvtxs  ~= length(vstat))
  msg = sprintf('nvtxs = %d, != length(vstat) = %d\n',nvtxs,length(vstat));
  err = 1;
  qoe(msg);error(msg);
end

tmp = [vtxno xyz vstat];

fid = fopen(filename,'w');
if(fid == -1)
  msg = sprintf('Could not open %s\n',filename);
  err = 1;
  qoe(msg);error(msg);
end

fprintf(fid,'%s\n',hdr);
fprintf(fid,'%d\n',nvtxs);
fprintf(fid,'%4d   %8.3f  %8.3f  %8.3f     %f\n',tmp'); %'

fclose(fid);

return;
