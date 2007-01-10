function [vtxno, xyz, vstat, hdr, msg] = fast_ldlabel(labelfile)
%
% [vtxno xyz vstat hdr msg] = fast_ldlabel(labelfile)
%
%


%
% fast_ldlabel.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
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

vtxno = [];
xyz = [];
vstat = [];
msg = 'OK';

if(nargin ~= 1) 
  msg = 'USAGE: [vtxno xyz vstat hdr msg] = fast_ldlabel(labelfile)';
  qoe(msg); error(msg);
end

labelfile = deblank(labelfile);
fid = fopen(labelfile);
if(fid == -1)
  msg = sprintf('Could not open  %s',labelfile);
  return;
end

hdr = fgetl(fid); % First line is like a header.
fgetl(fid);       % this is the number of vertices, redundant

m = fscanf(fid,'%f');

if(rem(length(m),5) ~= 0)
  msg = sprintf('Format of %s is bad',labelfile);
  return;
end

n = length(m)/5;

q = reshape(m, [5 n])'; %'

vtxno = q(:,1);
xyz   = q(:,[2:4]);
vstat = q(:,5);

fclose(fid);

return;
