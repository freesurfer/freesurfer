function [vtxno, xyz, vstat, hdr, msg] = fast_ldlabel(labelfile)
%
% [vtxno xyz vstat hdr msg] = fast_ldlabel(labelfile)
%
% $Id: fast_ldlabel.m,v 1.1 2003/03/04 20:47:38 greve Exp $

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
