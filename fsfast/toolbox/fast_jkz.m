function [jkz, p, fmn, fstd] = fast_jkz(f)
% [jkz, p, fmn, fstd] = fast_jkz(f)
%
% Computes the jackknifed z-score.
%
% jkz is the z-score of a column of f when the std is 
% computed by leaving that column out. 
%
% $Id: fast_jkz.m,v 1.1 2004/05/13 02:49:48 greve Exp $
%
% (c) Douglas N. Greve, 2004
%

if(nargin ~= 1)
  fprintf('[jkz, fmn, fstd, p] = fast_jkz(f)\n');
  return;
end

[ntp ns] = size(f);
fmn  = zeros(size(f));
fstd = zeros(size(f));

for jkslice = 1:ns
  ind = find([1:ns] ~= jkslice);
  fmn(:,jkslice)  = mean(f(:,ind),2);
  fstd(:,jkslice) = std(f(:,ind),[],2);
end 

jkz = (f-fmn)./fstd;

if(nargout < 2) return; end

p = erfc(jkz);

return;

