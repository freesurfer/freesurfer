function [ind1, ind2, w1, w2] = fidinterp(ped,toffset,echospacing,nechoes)
% [ind1, ind2, w1, w2] = fidinterp(ped,toffset,echospacing,nechoes)
%
% ped = post-excitation delay
% toffset = time to first FID echo
% echospacing = time between FID echoes
% nechoes = number of FID echoes
%
% ind1, ind2 = echo index
% w1, w2 = weights to apply to echo index
%
% fid(ped) = w1*fid(ind1) + w2*fid(2)
%
%


%
% fidinterp.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
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

ind1=[];
ind2=[];
w1=[];
w2=[];

if(nargin ~= 4)
  fprintf(' [ind1, ind2, w1, w2] = fidinterp(ped,toffset,echospacing,nechoes)\n');
  return;
end

ind = (ped-toffset)/echospacing + 1;
ind1 = floor(ind);
ind2 = ceil(ind);

itmp = find(ind1 < 1);
ind1(itmp) = 1;
itmp = find(ind1 > nechoes);
ind1(itmp) = nechoes;

itmp = find(ind2 < 1);
ind2(itmp) = 1;
itmp = find(ind2 > nechoes);
ind2(itmp) = nechoes;

w2 = ind-ind1;
w1 = 1-w2;

return;

