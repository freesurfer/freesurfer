function ypost = fmri_spatfil(ypre, SpatialFil)
%
% ypost = fmri_spatfil(ypre, SpatialFilter)
%
% Spatially filters a functional slice.
%
% -------- Input Arguments ------
% 1. ypre - functional slice to be processed. Dimension is
%       [nRows nCols nTP nRuns].
% 2. SpatialFilter - coefficients of spatial filter (2D).
%
% -------- Output Arguments ------
% 1. ypost - spatially filtererd slice of the same dimension
%       as the preprocessed slice, ie, [nRows nCols nTP nRuns].
%
% Douglas Greve (greve@nmr.mgh.harvard.edu)
% January 25, 1999
% February 12, 1999
%
%
%


%
% fmri_spatfilter.m
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

% Check that the number of input arguments is correct %
if(nargin ~= 2)
  msg = 'ypost = fmri_spatfil(ypre, SpatialFilter)';
  qoe(msg);error(msg);
end

% Get function slice dimensions %
nRows = size(ypre,1);
nCols = size(ypre,2);
nTP   = size(ypre,3);
nRuns = size(ypre,4);
nVoxels = nRows*nCols;

y = ones(nRows,nCols);
EdgeCorrection = conv2(y,SpatialFil,'same');

%%% Initialize matricies %%
ypost = zeros(nRows,nCols,nTP,nRuns);

%% Pass through each run %%
for r = 1:nRuns,

  %% Put slice in a temporary variable %%
  y = ypre(:,:,:,r);

  %% Pass through each time point %%
  for n = 1:nTP,
     z = conv2(y(:,:,n),SpatialFil,'same') ./ EdgeCorrection;
     y(:,:,n) = z;
   end

  ypost(:,:,:,r) = reshape(y, [nRows nCols nTP]);

end


return;
