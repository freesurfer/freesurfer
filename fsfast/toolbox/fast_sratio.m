function sratio = fast_sratio(num,den)
% sratio = fast_sratio(num,den)
% Computes a "signed" ratio, ie, 
%  if num > den then sratio = +num/den
%  if num < den then sratio = -den/num
%
% $Id: fast_sratio.m,v 1.2 2008/09/15 21:31:40 greve Exp $

%
% fast_fast_sratio.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2008/09/15 21:31:40 $
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

if(nargin ~= 2)
  fprintf('sratio = fast_sratio(num,den)\n');
  return;
end

mri = [];
if(isfield(num,'vol'))  
  num = num.vol; 
  mri = num;
end
if(isfield(den,'vol'))  
  mri = den;
  den = den.vol; 
end
 
sratio = zeros(size(num));
ind = find(num > den & den ~= 0);
sratio(ind) = num(ind)./den(ind);
ind = find(num < den & num ~= 0);
sratio(ind) = -den(ind)./num(ind);

if(~isempty(mri))
  mri.vol = sratio;
  sratio = mri;
end


return;


