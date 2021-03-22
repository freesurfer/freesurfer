function vvar = fmri_voxvar(e,dof,Cn,TPExcl)
%
% Computes the variance at a voxel given the error
% waveform at that voxel.
%
% vvar = fmri_voxvar(e,dof,<<Cn>,TPExcl>)
%
% e:   nRows x nCols x nTP x nRuns
% dof: scalar-int
% Cn:  nTPxnTP(xnRuns)
% TPExcl: nTPxnRuns
%
%


%
% fmri_voxvar.m
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

if(nargin ~= 2 & nargin ~= 3 & nargin ~= 4)
  msg = 'USAGE: fmri_voxvar(e,dof,<<Cn>,TPExcl>)';
  qoe(msg);
  error(msg);
end

if(nargin == 2)
  Cn = [];
  TPExcl = [];
end
if(nargin == 3)
  TPExcl = [];
end

sze = size(e);
ns = length(sze);

nRows = size(e,1);
nCols = size(e,2);
nTP   = size(e,3);
nRuns = size(e,4);

% Compute number of voxels %
nV = nRows*nCols; 

e = reshape(e, [nV nTP nRuns]);
e = permute(e, [2 1 3]);

vvar = 0;

if(isempty(Cn)) % Assume white noise, ie Cn = I
  for r = 1:nRuns,
     if(isempty(TPExcl))
       vvar = vvar + sum( e(:,:,r).^2 , 1);
     else
       ind = find(TPExcl(:,r) == 0);
       vvar = vvar + sum( e(ind,:,r).^2 , 1);
     end
  end

else % Use Cn as noise covariance matrix

  if(size(Cn,1) ~= nTP | size(Cn,2) ~= nTP | ...
     (size(Cn,3) ~=1 & size(Cn,3) ~= nRuns))
    msg = 'e and Cn dimensions are inconsistent';
    qoe(msg);
    error(msg);
  end

  if(size(Cn,3) == 1) % only one matrix specified, use for all runs
    iCn = inv(Cn);
    % chiCn = chol(iCn);
    for r = 1:nRuns,
      %ew = chiCn * e(:,:,r);
      ew =  e(:,:,r);
      vvar = vvar + sum(ew .* (iCn * ew),1);
    end

  else  % use different Cn for each run
    for r = 1:nRuns,
      iCn = inv(Cn(:,:,r));
      %ew = chol(iCn) * e(:,:,r);
      ew = e(:,:,r);
      vvar = vvar + sum(ew .* (iCn * ew),1);
      %vvar = vvar + ( sum(e(:,:,r) .* (inv(Cn(:,:,r)) * e(:,:,r)),1) );
    end  
  end

end

vvar = vvar/dof;


vvar = squeeze(reshape(vvar, [nRows nCols]));

return;
