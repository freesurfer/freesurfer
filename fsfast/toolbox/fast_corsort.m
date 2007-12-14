function [csort rhosort zrhoabs] = fast_corsort(rho,Nps)
% [csort rhosort zrhoabs] = fast_corsort(rho,Nps)
% rho is an NxN correlation matrix
% Nps is a list of number of inputs for each subject
%   sum(Nps) must equal N (the size of rho)
%   Eg, if there were 3 subjects each with 4 inputs, then
%      Nps = [4 4 4];
%   Eg, if the 2nd subject had 5 then
%      Nps = [4 5 4];
% csort will be Nmin by Ns, where Nmin = min(Nps) and Ns
%   is the number of subjects. Each row of csort represents
%   a "coherent component". The number in csort is the
%   subject component that belongs to that coherent component.
% rhosort - correlation values corresponding to the components
% zrhoabs - used for debugging.
%
% Eg:
%   y is component-by-space matrix
%   yn = fast_fnorm(y,2,1); % Normalize across space
%   rho = yn*yn'; % Cor Coeff across component
%
% $Id: fast_corsort.m,v 1.1 2007/12/14 19:19:37 greve Exp $

%
% fast_corsort.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/12/14 19:19:37 $
%    $Revision: 1.1 $
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

csort   = [];
rhosort = [];
zrhoabs = [];

if(nargin ~= 2)
  fprintf('[csort rhosort zrhoabs] = fast_corsort(rho,Nps)\n');
  return;
end

% Sort over absolute
rhoabs = abs(rho);

% Number of subjects
Ns = length(Nps);

% Tot number of comp across all subj
Ntot = sum(Nps); 
if(size(rho,1) ~= Ntot)
  fprintf('ERROR: dimension mismatch %d %d\n',size(rho,1), Ntot);
  return;
end

% Max number of components to sort
Nmax = min(Nps); 

% Create a list of subjects start and stop components
i1 = cumsum([0 Nps(1:end-1)]) + 1;
i2 = i1 + Nps - 1;

% Create lists for easy access
slist = zeros(Ntot,1); % subject at given row or col in rho
clist = zeros(Ntot,1); % component no at given row or col in rho
for s = 1:Ns
  slist(i1(s):i2(s)) = s;
  clist(i1(s):i2(s)) = 1:Nps(s);
end

% Zero values that are irrelevant. These are the values that
% are correlations within a given subject.
blk = [];
for s = 1:Ns
  bs = ones(Nps(s));
  blk = blkdiag(blk,bs);
end
aa = zeros(1,Ntot);
aa(1) = 1;
tt = toeplitz(ones(Ntot,1),aa);
indz = find(blk==1 | tt==1);
rhoabs(indz) = 0;
zrhoabs = rhoabs; % keep a copy

% Start loop over sorted components
csort  = zeros(Nmax,Ns);
rhosort = zeros(Nmax,Ns);
rhoabsR = rhoabs;
for nthcomp = 1:Nmax
  rhoabs = rhoabsR;

  %fprintf('%d ------------------------------------\n',nthcomp);
  %rhoabs

  while(1)
    [rhomax indmax] = max(rhoabs(:));
    [imax jmax] = ind2sub(size(rhoabs),indmax);
    sA = slist(imax);
    cA = imax - i1(sA) + 1;
    sB = slist(jmax);
    cB = jmax - i1(sB) + 1;
    if(csort(nthcomp,sA) == 0)
      csort(nthcomp,sA) = cA;
      rhosort(nthcomp,sA) = rho(imax,jmax);
    end
    if(csort(nthcomp,sB) == 0)
      csort(nthcomp,sB) = cB;
      rhosort(nthcomp,sB) = rho(imax,jmax);
    end

    indA = find(slist == sA);
    indB = find(slist == sB);
    z = zeros(size(rhoabs));
    z(indA,indB) = 1;
    ind = find(z);
    rhoabs(ind) = 0;
    
    indA = find(slist == sA & clist ~= cA);
    rhoabs(:,indA) = 0;
    rhoabs(indA,:) = 0;
    
    indB = find(slist == sB & clist ~= cB);
    rhoabs(:,indB) = 0;
    rhoabs(indB,:) = 0;

    if(length(find(rhoabs ~= 0)) == 0) break; end
  
  end

  for s = 1:Ns
    ind = find(slist == s & clist == csort(nthcomp,s));
    rhoabsR(:,ind) = 0;
    rhoabsR(ind,:) = 0;
  end
  
end



