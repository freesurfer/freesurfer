function [gcomp, theta, thetahat, Fsig, beta, Rrow] = tdr_ghostcomp(kpcn,epipar)
% [gcomp theta thetahat Fsig] = tdr_ghostcomp(kpcn,epipar)
%
% kpcn is (3,ncols,n1,n2,n3) with the rows NOT LR reversed
% gcomp is (ncols n1 n2 n3) 
% revenfixed = reven .* gcomp (careful with transpose)
%
% $Id: tdr_ghostcomp.m,v 1.4 2007/05/11 03:53:32 greve Exp $

order = 1;   % poly order
rthresh = 4; % relative threshold

gcomp = [];
if(nargin ~= 2)
  fprintf('[gkcomp theta thetahat Fsig] = tdr_ghostcomp(kpcn,epipar)\n');
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
%Rrow = transpose(inv(Frow));
Rrow = inv(Frow);

% Flip even rows left-right
kpcn(2,:,:,:) = flipdim(kpcn(2,:,:,:),2);

% kpcn is (3,ncols,n1,n2,n3) with the rows NOT LR reversed
kpcn = reshape(kpcn,[3 ncols n1*n2*n3]);
kpcn = permute(kpcn,[2 3 1]);

pcn1 = Rrow*kpcn(:,:,1);
pcn2 = Rrow*kpcn(:,:,2);
pcn3 = Rrow*kpcn(:,:,3);

kpcnref = (kpcn(:,:,1) + kpcn(:,:,3))/2;
kpcnmov = kpcn(:,:,2);

kpcnaref = (abs(kpcn(:,:,1)) + abs(kpcn(:,:,3)))/2;
kpcnamov = abs(kpcn(:,:,2));

% Reconstruct the rows
%pcnref = (Rrow*kpcnref);
%pcnmov = (Rrow*kpcnmov);
pcnref = (pcn1+pcn3)/2;
pcnmov = pcn2;


[thetahat theta Fsig beta] = tdr_phdist_est(pcnref,pcnmov,rthresh,order);

gcomp = exp(+i*thetahat);

% Apply compensation to the pcnmov
pcnmovcomp = pcnmov .* gcomp;

% Convert it back to kspace
kpcnmovcomp = Frow*pcnmovcomp;

% The compensated should look very close to the ref
nn = 1:ncols;
plot(nn,kpcnaref(:,1),'+-',nn,kpcnamov(:,1),nn,abs(kpcnmovcomp(:,1)));
legend('ref','mov','movcomp')

gcomp = reshape(gcomp,[ncols n1 n2 n3]);

gcompmov = exp(+i*thetahat/2);
gcompref = exp(-i*thetahat/2);

return;

%
% tdr_ghostcomp.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2007/05/11 03:53:32 $
%    $Revision: 1.4 $
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

