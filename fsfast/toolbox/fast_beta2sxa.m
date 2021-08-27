function [hsxa, hCovMtx] = fast_beta2sxa(beta,rvar,Nc,Npc,X)
% [hsxa, hCovMtx] = fast_beta2sxa(beta,rvar,Nc,Npc,X)
% Converts beta/rvar into selxavg format
% beta full reg coeff set
% rvar resid error variance
% Nc - number of non-null conditions
% Npc - number of reg per cond
% X - design matrix (multiplied by W if whiten)


%
% fast_beta2sxa.m
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

if(nargin ~= 5)
  fprintf('[hsxa, hCovMtx] = fast_beta2sxa(beta,rvar,Nc,Npc,X)\n');
  return;
end

nTask = Nc*Npc;	
nv = size(beta,2);
XtX0 = X'*X;
iXtX = inv(XtX0);
hCovMtx = iXtX(1:nTask,1:nTask);
XtX = XtX0(1:nTask,1:nTask);

hhattmp = beta(1:nTask,:); % Remove nuisance
hhattmp = [zeros(Npc,nv); hhattmp]; % Add zero for cond 0
hhattmp2 = reshape(hhattmp,[Npc Nc+1 nv]);

hstd = sqrt( (diag(hCovMtx).*diag(XtX)) * rvar);
hstdtmp = hstd(1:nTask,:); % Remove offset and baseline
hstdtmp = [repmat(sqrt(rvar), [Npc 1]); hstdtmp]; % Add 0 for cond 0
hstdtmp2 = reshape(hstdtmp,[Npc Nc+1 nv]);

%--- Merge Averages and StdDevs ---%
hsxa = zeros(Npc,2,Nc+1,nv);
hsxa(:,1,:,:) = hhattmp2;
hsxa(:,2,:,:) = hstdtmp2;

return;
