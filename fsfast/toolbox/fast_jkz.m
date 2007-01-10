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

