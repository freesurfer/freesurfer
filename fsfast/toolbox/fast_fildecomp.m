function [U, S, F, pve0, niters, err, TO] = ...
      fast_fildecomp(Forig,cpvemin,errtol,nitersmax)
%
% [U S F pve0 niters err TO] = fast_fildecomp(Forig,cpvemin,errtol,nitersmax)
%
% 


%
% fast_fildecomp.m
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

U = [];
S = [];
F = [];
niters = 0;
err = 10^10;
TO = 0;

if(nargin ~= 4)
  fprintf('[U S F niters err TO] = fast_fildecomp(Forig,cpvemin,errtol,nitersmax)\n');
  return;
end

F = Forig;
[Uorig Sorig tmp] = svd(Forig);
N = size(F,1);
c = ones(N,1);

while(err > errtol & niters <= nitersmax)

  F = F ./ repmat(c,[1 N]);
  [U S tmp] = svd(F);
  ds = diag(S);
  cpve = 100*(cumsum(ds)/sum(ds));
  nkeep = min(find(cpve>=cpvemin));
  nn = 1:nkeep;
  U = U(:,nn);
  S = S(nn,nn);
  F = U*S*U'; %'
  c = sum(F,2);
  err = std(c);
  niters = niters + 1;

end

if(niters == nitersmax) TO = 1; end

Snew = U'*Forig*U; %'
pve0 = 100*sum(diag(Snew))/sum(diag(Sorig));

return;
