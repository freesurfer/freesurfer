function [pSigMax, vpSigMax] = fmri_pmax(pSig,vSig)
%
% [pSigMax vpSigMax] = fmri_pmax(pSig)
% [pSigMax vpSigMax] = fmri_pmax(pSig, vSig)
%
% Returns the maximum significance from a vector of significance
% values.  pSig may be a 3d matrix in which case the maximum is
% searched for over the third dimension.  pSigMax would then be
% a 2d matrix of maximum significances at each point, corrected
% for the number of unplanned comparisons (Bonferoni correction).
%
% If vSig is included, it must be of the same dimension as pSig.
% vpSigMax is a 2d matrix with values from vSig that correspond
% to the same locations in pSig from which the values in pSigMax 
% came.
%
%


%
% fmri_pmax.m
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

if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE: [pSigMax <vpSigMax>] = fmri_pmax(pSig, <vSig>)';
  qoe(msg);error(msg);
end

if(nargout == 2 & nargin == 1)
  msg = 'Must specify vSig to get vpSigMax';
  qoe(msg);error(msg);
end

%% Get dimensions %%
[r c nt] = size(pSig);
nv = r*c; % number of voxels %

%% reshape into time-points by voxels %%
pSig = reshape(pSig, [nv nt])'; %'

%% Get the maximum p values for each voxel %%
[pSigMax nt_pmax] = max(pSig);

% Bonferoni Correction %
pSigMax = pSigMax/nt;

% Put back into rows and colums %
pSigMax = reshape(pSigMax,[r c]);

% Return if vpSigMax is not needed %
if(nargout == 1) return; end


% Get vSig values corresponding to pSigMax %
vSig = reshape(vSig, [nv nt])'; %'

%% Compute indicies of pSigMax voxel-timepoints %%
i_pmax = [0:nv-1]*nt + nt_pmax;

%% Extract vSig values %%
vpSigMax = vSig(i_pmax);

% Put them back into rows and colums %
vpSigMax = reshape(vpSigMax,[r c]);


return;
