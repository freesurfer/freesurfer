function w = fisher_twoclass(D1,D2)
% w = fisher_twoclass(D1,D2)
%
% Fisher two-class linear discriminant
%
% D1 = n1 by nfeatures
% D2 = n2 by nfeatures
%
% n1 = number of samples of group 1
% n2 = number of samples of group 2
% n1+n2 must be > nfeatures
%
% w is the nfeatures-by-1 optimal weighting
%
% To test:
%   D1 = randn(100000,2) + 1;
%   D2 = randn(100000,2) - 1;
%   w = fisher_twoclass(D1,D2)
%   w should be close to [0.707 0.707] (or neg of that)
%

%
% fisher_twoclass(D1,D2)
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

w=[];
if(nargin ~= 2)
  fprintf('w = fisher_twoclass(D1,D2)\n');
  return;
end

if(size(D1,2) ~= size(D2,2))
  fprintf('ERROR: dimension mismatch\n');
  return;
end

nfeatures = size(D1,2);

n1 = size(D1,1);
n2 = size(D2,1);
ntot = n1 + n2;

if(nfeatures >= ntot) 
  fprintf('ERROR: nfeatures (%d) >= ntot (%d)\n',nfeatures,ntot);
  return;
end

D1mn = mean(D1); % mean feature vector, group 1
R1 = D1 - repmat(D1mn,[n1 1]);
S1 = R1'*R1;

D2mn = mean(D2); % mean feature vector, group 2
R2 = D2 - repmat(D2mn,[n2 1]);
S2 = R2'*R2;

Sw = S1 + S2;

DmnDiff = D1mn - D2mn;
Sb = DmnDiff'*DmnDiff;

M = inv(Sw)*Sb;
[evects evals] = eig(M);
[emax imax] = max(diag(evals));
w = evects(:,imax);

return;

