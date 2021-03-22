function [jkz, p, fmn, fstd] = fast_jkz(f)
% [jkz, p, fmn, fstd] = fast_jkz(f)
%
% Computes the jackknifed z-score.
%
% jkz is the z-score of a column of f when the std is 
% computed by leaving that column out. 
%
%
%
% (c) Douglas N. Greve, 2004
%


%
% fast_jkz.m
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

