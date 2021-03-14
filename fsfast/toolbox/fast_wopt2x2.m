function [w, lambda] = fast_wopt2x2(Cs,Cr)
% [w lambda] = fast_wopt2x2(Cs,Cr)
%
% Cs is the 3xN Signal Covariance Matrix
%
% Cs(n) = [Cs(1,n) Cs(3,n);
%          Cs(3,n) Cs(2,n)];
%
% Cs(1,:) = sum(s1.^2);
% Cs(2,:) = sum(s2.^2);
% Cs(3,:) = sum(s1.*s2);
%
% Cr is the 3xN Noise Covariance Matrix
%
% Cr(n) = [Cr(1,n) Cr(3,n);
%          Cr(3,n) Cr(2,n)];
%
% Cr(1,:) = sum(r1.^2);
% Cr(2,:) = sum(r2.^2);
% Cr(3,:) = sum(r1.*r2);
% 
% lambda is the 1xN maximum eigval of inv(Cr)*Cs
% w is the corresponding 2xN eigvect
%


%
% fast_wopt2x2.m
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

w = [];

if(nargin ~= 2)
  fprintf('[w lambda] = fast_wopt2x2(Cs,Cr)\n');
  return;
end

nv = size(Cs,2);
if(size(Cr,2) ~= nv)
  fprintf('ERROR: dimension mismatch\n');
  return;
end

% Determinant of Cr %
detCr = Cr(1,:).*Cr(2,:)-Cr(3,:).^2;

% Check that all Cr are postive def %
indneg = find(detCr <= 0);
nindneg = length(indneg);
if(nindneg ~= 0)
  fprintf('ERROR: found %d non-positive definite matrices.\n', ...
	  nindneg);
  fprintf('       First at index %d\n',indneg(1));
end

% Compute inv(Cr)*Cs = [a b; c d]; %
% Note: in general, not a sym matrix %
CrCs3 = Cr(3,:).*Cs(3,:);
a = (Cr(2,:).*Cs(1,:) - CrCs3)./detCr;
b = (Cr(2,:).*Cs(3,:) - Cr(3,:).*Cs(2,:))./detCr;
c = (Cr(1,:).*Cs(3,:) - Cr(3,:).*Cs(1,:))./detCr;
d = (Cr(1,:).*Cs(2,:) - CrCs3)./detCr;

% Maximum eigenvalue %
lambda = ((a+d) + sqrt( (a+d).^2 + 4*(b.*c - a.*d)))/2;

% Corresponding eigenvector %
w(1,:) = ones(1,nv);
w(2,:) = -(a-lambda)./b;
w = w./repmat(sqrt(sum(w.^2)),[2 1]);

return;

















