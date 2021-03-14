function [Nc, CondList, Holes] = fast_par2nconds(par)
% [Nc, CondList, Holes] = fast_par2nconds(par)
%
% Returns the number of conditions, excluding 0. CondList
% is the list of conditions found. If the condition numbers
% are not contiguous starting at 1, then Holes = 1. If there
% are holes, then returns max(CondList) instead of 
% length(CondList). If there are no holes, then these are the
% same.
%
% par is (Nstim,2,Nruns) 


%
% fast_par2nconds.m
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

Nc = [];
CondList = [];
Holes = [];

if(nargin ~= 1)
  fprintf('[Nc, CondList, Holes] = fast_par2nconds(par)\n');
  return;
end

cids = par(:,2,:);
inz = find(cids ~= 0);
cids = cids(inz);

CondList = unique(cids);
Nc = max(CondList);

d = diff(CondList);
nn1 = length(find(d ~= 1));
if(nn1 ~= 0) Holes = 1;
else         Holes = 0;
end

if( CondList(1) ~= 1 ) Holes = 1; end

return;
