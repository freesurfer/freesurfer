function y = fmri_shiftcol(x,n,wrapflag)
%
% y = fmri_shiftcol(x,n)
% y = fmri_shiftcol(x,n,wrapflag)
%
% Shift the columns of x by n.  If n > 0, then
% the shift is to the right (ie, the first column
% of x becomes the (n+1)th column of y, the second
% column of x becomes the (n+2)th column of y, etc).
% If n < 0, the columns are shifted to the left.
% If wrapflag = 1 or is not specified, the columns
% are wrapped around, otherwise they are filled with
% zeros.  x must be a 2D array.  To shift rows, pass
% the transpose of x.
%
% $Id: fmri_shiftcol.m,v 1.1 2003/03/04 20:47:40 greve Exp $

if(nargin ~= 2 & nargin ~= 3)
  msg = 'Usage: y = fmri_shiftcol(x,n,<wrapflag>)';
  qoe(msg); error(msg);
end

if(length(size(x)) ~= 2)
  msg = 'x must be a 2D array';
  qoe(msg); error(msg);
end

if(nargin == 2) wrapflag = 1; end

if(n==0) 
  y = x;
  return;
end

ncx = size(x,2);

if(wrapflag ~= 0)  x = [x x];
else 
  if(n>0)  x = [zeros(size(x)) x];
  else     x = [x zeros(size(x))];
  end
end

n = mod(-n+ncx,ncx);
y = x(:,n+1:n+ncx);

return;
