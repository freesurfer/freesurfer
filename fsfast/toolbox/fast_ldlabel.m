function [vtxno, xyz, vstat, hdr, msg] = fast_ldlabel(labelfile)
%
% [vtxno xyz vstat hdr msg] = fast_ldlabel(labelfile)
%
%


%
% fast_ldlabel.m
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
