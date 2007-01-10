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
% '$Id: fmri_mrestriction2.m,v 1.2 2007/01/10 22:02:33 nicks Exp $'


%
% fmri_mrestriction2.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
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
