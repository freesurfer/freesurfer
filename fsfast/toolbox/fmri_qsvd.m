function [Uhat, Shat, Sall] = fmri_qsvd(D,cutofftype,cutoffthresh)
%
% [Uhat, Shat, Sall] = fmri_qsvd(D,cutofftype,cutoffthresh)
%
% D - 2d matrix
% cutofftype, cutoffthreshold
%   lmin   - exclude eigenvalues below the cutoffthreshold
%   nmax   - include only cutoffthreshold eigenvalues
%   pctvar - include enough eigenvalues to explain cutoffthreshold
%            percent of the variance.  Range: 0-1.
%
% Uhat - Matrix whose columns are eigenvectors.
% Shat - vector of corresponding singular/eigen values 
% Sall - vector of all singular/eigen values 
%
% Computes [U S2 V] = svd(D'*D), and keeps the appropriate number
% of eigencomponents.  Note that U = V and that S2 are the square
% of the eivenvalues. sqrt(S2) is returned.
%
%


%
% fmri_qsvd.m
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

if(nargin ~= 3)
  msg = 'USAGE: [Uhat, Shat] = fmri_qsvd(D,cutofftype,cutoffthresh)';
  qoe(msg); error(msg);
end

cutofftype = lower(cutofftype);
if(~strcmp(cutofftype,'lmin') & ...
   ~strcmp(cutofftype,'nmax') & ...
   ~strcmp(cutofftype,'pctvar') )
  msg = sprintf('cutofftype = %s, must be lmin, nmax, or pctvar',cutofftype);
  qoe(msg); error(msg);
end

if(length(size(D)) ~= 2)
  msg = 'D must be a 2-D matrix';
  qoe(msg); error(msg);
end

Dmaxrank = min(size(D));

[U S tmp] = svd(D'*D/size(D,1));
Sall = diag(S);
clear tmp S;

switch(cutofftype)

  case {'lmin'},
    lmin = cutoffthresh;
    ind = find(Sall > lmin);

  case {'nmax'},
    nmax = cutoffthresh;
    m = min(nmax,Dmaxrank);
    ind = [1:m];

  case {'pctvar'},
    pctvar = cutoffthresh;
    if(pctvar > 1) pctvar = pctvar/100; end
    m = min(find(cumsum(Sall)/sum(Sall) > pctvar));
    ind = [1:m];

end

Sall = sqrt(Sall);
Shat = Sall(ind);
Uhat = D*U(:,ind)*diag(1./Shat)/sqrt(size(D,1));

return;
