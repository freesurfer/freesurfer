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

