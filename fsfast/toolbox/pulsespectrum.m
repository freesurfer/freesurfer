function f = pulsespectrum(x0,x1,w)
% f = pulsespectrum(x0,x1,w)
%
% Continuous (analytical) spectrum of a pulse starting 
% at x0 and ending at x1. w is the radial frequency.
%
% If x0 and x1 have more than one element, then a spectrum
% is produced, and f will have multiple columns.
% 
% If x0 = -x1, then the result will be all real.
%
% $Id: pulsespectrum.m,v 1.1 2004/01/17 05:32:45 greve Exp $

if(nargin ~= 3)
  fprintf('f = pulsespectrum(x0,x1,w)\n');
  return;
end

if(length(x0) ~= length(x1))
  fprintf('ERROR: x0 and x1 have different sizes\n');
  return;
end

d = x1-x0;
ind = find(d <= 0);
if(~isempty(ind))
  fprintf('ERROR: x0 >= x1\n');
  return;
end
nx = length(x1);

indwz  = find(w==0);
indwnz = find(w~=0);

f = zeros(length(w),nx);
for nthx = 1:nx
  x00 = x0(nthx);
  x11 = x1(nthx);

  fx = zeros(size(w));
  fx(indwz) = x11-x00;

  wx0 = w(indwnz)*x00;
  wx1 = w(indwnz)*x11;
  fx(indwnz) = (cos(wx1) - cos(wx0) + i*(sin(wx0)-sin(wx1)) )./(-i*w(indwnz));

  f(:,nthx) = fx;
end

return
