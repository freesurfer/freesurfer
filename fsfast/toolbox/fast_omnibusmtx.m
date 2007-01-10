function R = fast_omnibusmtx(nEventTypes,TER,TW,TPS,RmPrestim)
% 
% R = fast_omnibusmtx(nEventTypes,TER,TW,TPS,RmPrestim)
%
% Allows computation of an omnibus matrix that includes
% zeroing the prestimulus baseline and computing and 
% computing the stat based on the poststim component.
%
% 


%
% fast_omnibusmtx.m
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

R = [];

if(nargin ~= 4 & nargin ~= 5)
  fprintf('USAGE: R = fast_omnibusmtx(nEventTypes,TER,TW,TPS,RmPrestim)');
  return;
end

if(nargin == 4) RmPrestim = 0; end

Nh = round(TW/TER);
NhTot = Nh * nEventTypes;

if(isempty(TPS)) TPS = 0; end

Nps = round(TPS/TER) + 1; % +1 inclues 0

if(RmPrestim)
  Rps = -ones(Nh-Nps,Nps)/Nps;
  %R1 = [Rps eye(Nh-Nps,Nh) zeros(Nh-Nps,(nEventTypes-1)*Nh)];
  R1 = [Rps eye(Nh-Nps,Nh-Nps)];
  R0 = zeros(size(R1));
  R = [];
  for n = 1:nEventTypes;
    Ra = repmat(R0,[1 (n-1)]);
    Rz = repmat(R0,[1 (nEventTypes-n)]);
    Rc = [Ra R1 Rz];
    R = [R; Rc];
  end
else
  R = eye(NhTot);
end

return;
