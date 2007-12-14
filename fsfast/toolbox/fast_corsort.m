function [csort rhosort spair zrhoabs] = fast_corsort(rho,Nps)
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
%   Each row will have all the subjects; each col must not have
%   any repeats.
% rhosort - correlation values corresponding to the components
% spair - Nmin by Ns. The value is the subject that the column
%   subject was paired with.
% zrhoabs - used for debugging.
%
% Eg:
%   y is component-by-space matrix
%   yn = fast_fnorm(y,2,1); % Normalize across space
%   rho = yn*yn'; % Cor Coeff across component
%
% $Id: fast_corsort.m,v 1.2 2007/12/14 20:11:19 greve Exp $

%
% fast_corsort.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/12/14 20:11:19 $
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

csort   = [];
rhosort = [];
zrhoabs = [];
spair   = [];

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

% Start loop over coherent components
csort  = zeros(Nmax,Ns);
rhosort = zeros(Nmax,Ns);
rhoabsR = rhoabs;
for nthcomp = 1:Nmax
  rhoabs = rhoabsR;

  % Loop over subjects
  nth = 0;
  while(1)
    nth = nth + 1;

    % Find max of the (remaining) cors
    [rhomax indmax] = max(rhoabs(:));

    % Extract the corresponding subject and components numbers
    [imax jmax] = ind2sub(size(rhoabs),indmax);
    sA = slist(imax);
    cA = imax - i1(sA) + 1;
    sB = slist(jmax);
    cB = jmax - i1(sB) + 1;

    % Add them to the list
    if(csort(nthcomp,sA) == 0)
      csort(nthcomp,sA) = cA;
      rhosort(nthcomp,sA) = rho(imax,jmax);
      spair(nthcomp,sA) = sB;
    end
    if(csort(nthcomp,sB) == 0)
      csort(nthcomp,sB) = cB;
      rhosort(nthcomp,sB) = rho(imax,jmax);
      spair(nthcomp,sB) = sA;
    end

    % Eliminate the sA-sB block 
    indA = find(slist == sA);
    indB = find(slist == sB);
    z = zeros(size(rhoabs));
    z(indA,indB) = 1;
    ind = find(z);
    rhoabs(ind) = 0;
    
    % For each subject, eliminate the components NOT
    % chosen for further consideration for this
    % coherent component
    indA = find(slist == sA & clist ~= cA);
    rhoabs(:,indA) = 0;
    rhoabs(indA,:) = 0;
    indB = find(slist == sB & clist ~= cB);
    rhoabs(:,indB) = 0;
    rhoabs(indB,:) = 0;

    % Now eliminate now eliminate cross blocks of
    % subjects that are neither sA or sB
    if(nth == 1)
      ind = find(slist ~= sA & slist ~= sB);
      rhoabs(ind,ind) = 0;
    end

    % Is there anything left?
    if(length(find(rhoabs ~= 0)) == 0) break; end
  
  end

  % Now eliminate the components just selected from
  % the master list
  for s = 1:Ns
    ind = find(slist == s & clist == csort(nthcomp,s));
    rhoabsR(:,ind) = 0;
    rhoabsR(ind,:) = 0;
  end
  
end



