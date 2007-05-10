function [gcomp, theta, thetahat, Fsig] = tdr_ghostcomp(kpcn,epipar)
% [gcomp theta thetahat Fsig] = tdr_ghostcomp(kpcn,epipar)
%
% kpcn is (3,ncols,n1,n2,n3) with the rows NOT LR reversed
% gcomp is (ncols n1 n2 n3)
% kevenfixed = keven .* gcomp (careful with transpose)
%
% $Id: tdr_ghostcomp.m,v 1.3 2007/05/10 23:55:16 greve Exp $

order = 1;   % poly order
rthresh = 4; % relative threshold

gcomp = [];
if(nargin ~= 2)
  fprintf('[gcomp theta thetahat Fsig] = tdr_ghostcomp(kpcn,epipar)\n');
  return;
end

[npcn ncols n1 n2 n3] = size(kpcn);
if(npcn ~= 3) 
  fprintf('ERROR: npcn = %d, should be 3)\n',npcn);
  return;
end

echospacing = epipar.echospacing;
tDwell      = epipar.tDwell;
tRampUp     = epipar.tRampUp;
tFlat       = epipar.tFlat;
tRampDown   = epipar.tRampDown;
tDelSamp    = epipar.tDelSamp;


[kvec gvec] = kspacevector2(ncols,tDwell,tRampUp,tFlat,...
			    tRampDown,tDelSamp,0);

% Compute the Ideal col and row DFT reconstruction matrices
Frow = fast_dftmtx(kvec);
Frow = fast_svdregpct(Frow,90);
Rrow = transpose(inv(Frow));


% Flip rows left-right
kpcn(2,:,:,:) = flipdim(kpcn(2,:,:,:),2);

% kpcn is (3,ncols,n1,n2,n3) with the rows NOT LR reversed
kpcn = reshape(kpcn,[3 ncols n1*n2*n3]);
kpcn = permute(kpcn,[2 3 1]);
kpcnref = (kpcn(:,:,1) + kpcn(:,:,2))/2;
kpcnmov = kpcn(:,:,2);

pcnref = (Rrow*kpcnref);
pcnmov = (Rrow*kpcnmov);

[thetahat theta Fsig beta] = tdr_phdist_est(pcnref,pcnmov,rthresh,order);

gcomp = exp(+i*thetahat);
gcomp = reshape(gcomp,[ncols n1 n2 n3]);


return;

%
% tdr_ghostcomp.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/05/10 23:55:16 $
%    $Revision: 1.3 $
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

