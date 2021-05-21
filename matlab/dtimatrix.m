function B = dtimatrix(bvalues,bvectors)
% B = dtimatrix(bvalues,bvectors)
%
% Constructs a DTI design matrix from the bvalues and bvectors.
% bvalues is a vector of length N
% bvectors is a matrix either Nx3 or 3xN
%
% B will be N by 7
% The 7th is the mean (all ones)
% The tensor will be constructed using the following regressors
%     1 2 3
%     2 4 5
%     3 5 6
%

%
% dtimatrix.m
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


if(nargin ~= 2)
  fprintf('B = dtimatrix(bvalues,bvectors)\n');
  return;
end

if(size(bvectors,2) ~= 3) bvectors = bvectors'; end
if(size(bvectors,2) ~= 3) 
  fprintf('ERROR: bvectors must be Nx3 or 3xN\n');
  return;
end
Nb = size(bvectors,1);

if(Nb ~= length(bvalues))
  fprintf('ERROR: dimension mismatch between bvectors and bvalues\n');
  return;
end

bvalues = bvalues(:);

B = zeros(Nb,7);
B(:,1) =     bvalues .* bvectors(:,1).*bvectors(:,1);
B(:,2) = 2 * bvalues .* bvectors(:,1).*bvectors(:,2);
B(:,3) = 2 * bvalues .* bvectors(:,1).*bvectors(:,3);
B(:,4) =     bvalues .* bvectors(:,2).*bvectors(:,2);
B(:,5) = 2 * bvalues .* bvectors(:,2).*bvectors(:,3);
B(:,6) =     bvalues .* bvectors(:,3).*bvectors(:,3);
B(:,7) = 1;

return;
