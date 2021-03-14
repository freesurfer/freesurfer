function [ySA, dofSA] = fmri_sxa2sa(eVar,Ch,hAvg,hd)
%
% [ySA dofSA] = fmri_sxa2sa(eVar,Ch,hAvg,hd)
%
% Converts from selxavg format to selavg format.  In selavg
% format, the averages and standard deviations are interleaved
% in the same data structure for all conditions, including
% fixation.  In selxavg, the average for the fixation condition
% is always zero (by hypothesis) with covariance matrix equal
% to eVar*I (ie, the variances are equal across delays, and
% the delays are not correlated with each other).
%
% Note: the information about Ch is lost in that it cannot
% be recovered from ySA and/or dofSA, it is saved elsewhere.
%
%


%
% fmri_sxa2sa.m
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


if(nargin ~= 4)
  msg = 'Incorrect number of input arguments';
  qoe(msg);
  error(msg);
end

[hd.Nrows hd.Ncols] = size(eVar);

Nv  = hd.Nrows*hd.Ncols;
Nch = hd.Nnnc * hd.Nh;

% DOFs %
vdof = diag(hd.SumXtX);

% Compute the standard deviation of each component %
hStd = sqrt( diag(Ch).*vdof * reshape1d(eVar)'); %'

% Reshape into a more convenient form %
hAvg = reshape(hAvg, [Nv Nch])'; %'
hAvg = reshape(hAvg, [hd.Nh hd.Nnnc hd.Nrows hd.Ncols]);
hStd = reshape(hStd, [hd.Nh hd.Nnnc hd.Nrows hd.Ncols]);

% Allocate space of selavg structure %
ySA = zeros(hd.Nh,2,hd.Nc,hd.Nrows,hd.Ncols);

% Condition 0: mean=0, std=eStd %
ySA(:,2,1,:,:) = repmat(reshape(sqrt(eVar),[1 hd.Nrows hd.Ncols]),[hd.Nh 1 1]);

% Fill in the avg and std for the other conditions %
ySA(:,1,[2:hd.Nnnc+1],:,:) = hAvg;
ySA(:,2,[2:hd.Nnnc+1],:,:) = hStd;

% ySA now has dimension [Nh 2 Nc Nrows Ncols]
% Reshape into the [Nrows Ncols Nh*2*Nc ]
nSAT = 2 * hd.Nc * hd.Nh;
ySA = reshape(ySA, [nSAT hd.Nrows hd.Ncols]); 
ySA = permute(ySA, [2 3 1]);

% Compute the degrees of freedom for each condition %
if(hd.GammaFit ~= 1)
  dof = diag(hd.SumXtX)'; %'
  dof = dof(1:hd.Nh:Nch);
  dof0 = hd.Ntp*hd.Nruns - sum(dof);
  dof = [dof0 dof];
else
  dof = hd.Npercond;
end
dofSA = dof;

return;
