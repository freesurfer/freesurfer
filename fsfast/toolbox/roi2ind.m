function ind = roi2ind(siz,roi)
%
% ind = roi2ind(siz,roi)
%
% Converts a Region of Interest into indicies. 
%
%   siz = [nrows ncols]
%   roi = [rowmin colmin rowmax colmax]
%
% $Id: roi2ind.m,v 1.1 2003/03/04 20:47:41 greve Exp $

if(nargin ~= 2)
  msg = 'USAGE: ind = roi2ind(siz,roi)';
  qoe(msg);error(msg);
end

r1 = roi(1);
c1 = roi(2);
r2 = roi(3);
c2 = roi(4);
nr = siz(1);
nc = siz(2);

if(r1 > r2) 
  tmp = r2;
  r2 =  r1;
  r1 = tmp;
end
if(c1 > c2) 
  tmp = c2;
  c2 =  c1;
  c1 = tmp;
end

if(nr <= 0 | r1 <= 0 | r2 > nr | nc <= 0 | c1 <= 0 | c2 > nc)
  msg = sprintf('Invalid Range: siz = %d %d, roi = %d %d %d %d',siz,roi);
  qoe(msg);error(msg);
end

rl = [r1:r2]';
cl = [c1:c2];

rll = repmat(rl, [1 length(cl)]);
cll = repmat(cl, [length(rl) 1]);

ind = reshape1d(((cll-1)*nr + rll)');  % 1-Based Column Major !

return;