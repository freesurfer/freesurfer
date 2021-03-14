function RM = fmri_mrestriction2(wcond,wdelay,sumconds,sumdelays)
% RM = fmri_mrestriction2(wcond,wdelay,sumconds,sumdelays)
%
% Creates a restriction matrix given the vector of condition 
% weights and the vector of delays weights.  
%
% Notes:
%  1. wcond is the vector of condition weightings NOT including
%     the null condition.
%  2. The number of rows in  RM will be:
%            1    for sumconds=1 and sumdelays=1
%       nconds    for sumdelays=1
%       ndelays   for sumconds=1
%       nc*nd     for sumconds=0 and sumdelays=0
%


%
% fmri_mrestriction2.m
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

if(nargin ~= 4)
  msg = 'USAGE: RM = fmri_mrestriction2(wcond,wdelay,sumconds,sumdelays)';
  qoe(msg);error(msg);
end

nconds  = length(wcond);  % number of NON-NULL conditions
ndelays = length(wdelay);
ncd     = ndelays*nconds;

% make sure they are row vectors %
wcond  = reshape(wcond,  [1 nconds]);
wdelay = reshape(wdelay, [1 ndelays]);

if(sumconds & sumdelays)
  RM = [];
  for c = 1:nconds
    RM = [RM wcond(c)*wdelay];
  end
  return;
end

if(sumconds)
  RM = [];
  n = 1;
  for d = 1:ndelays
    if(wdelay(d) ~= 0) 
      %fprintf('d= %2d , wdelay(d) = %g\n',d,wdelay(d));
      v = zeros(nconds,ndelays);
      v(:,d) = wdelay(d)*wcond'; %'
      v2 = reshape(v',[1 ncd]); %'
      RM(n,:) = v2;
      n = n + 1;
    end
  end
  return;
end

if(sumdelays)
  RM = [];
  v0 = zeros(nconds,ndelays);
  v0(1,:) = wdelay;
  for c = 1:nconds
    if(wcond(c) ~= 0) 
      v = wcond(c)*v0;
      v = fmri_shiftcol(v',c-1)';
      RM = [RM v];
    end
  end
  return;
end

% neither sumdelays nor sumconds %
RM = diag(reshape1d(wdelay'*wcond)); %'

return;
