function cdf = rft_zcluster_cdf(csize,zthresh,fwhm,ssize,D)
% cdf = rft_zcluster_cdf(csize,zthresh,fwhm,ssize,D)
%
% Prob that a cluster >= csize defined by threshold zthresh will be
% found in a D-dim z-field with fwhm smoothness of ssize search space.
% zthresh is the voxel-wise z thresh. ssize and csize are measured
% in non-resel units. csize, fwhm, and ssize are measured in the
% same units. This has function has been tested with simulated data
% for both D=2 and D=3. For D=2, the test was to use a single slice
% of data, not a cortical surface; as far as I can tell, this does
% not work for surface data. KJW's SurfStat does work for cortical
% surfaces, but it is not clear what it is using exactly (maybe
% taking topology into account).
%
% Based on:
% Friston, Worsley, Frackowiak, Mazziotta, Evans. Assessing the
% significance of focal activations using their spatial extent. HBM
% 1994, 1:214-220. 
%
% This function generates the results as found in Table 1 if the
% special case for D=3 is not used. See the code. If the special case
% for D=3 is used (the default), then it matches results from KJW's
% stat_volume.m and FSL and fits simulations pretty well (note:
% when clustering real data, need 26 connectivity, ie, faces,
% edges, and corners). 
%
% Also in:
% Friston, Holmes, Poline, Price, and Frith. Detecting
% Activations in PET and fMRI: Levels of Inference and Power.
% Neuroimage 40, 223-235 (1996).
% 

%
% rft_zcluster_cdf.m
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

cdf = [];
if(nargin ~= 5)
  fprintf('cdf = rft_zcluster_cdf(csize,zthresh,fwhm,ssize,D)\n');
  return;
end

% Variables names and equations are from Friston, et al.

% Height threshold meausred in z-units
u = zthresh; 

% Equivalent p-value threshold. Note that the paper uses phi(-u),
% but the results dont work out that way. Note that phi(-u) can be
% used if (1-phiu) is used instead of phiu in the equation for beta
% below. See Hayasaka and Nichols 2003, App A, equation 4.
phiu = fast_z2p(u);

k = csize; % cluster size to test (actual units, not resels)
S = ssize; % search space (actual units, not resels)

W = fwhm/sqrt(4*log(2));

% Expected number of clusters (Eq 2)
% This form appears to go back to Hasofer 1978
Em = exp(-(u.^2)/2) .* u.^(D-1) * (2*pi).^(-(D+1)/2) .* S ./ (W.^D);
if(D==3)
  % The only difference here is (u.^(D-1)-1), which appears to come
  % from Worsely 1996. This is the FSL implementation. This also
  % matches results from KJW's stat_volume.m
  Em = exp(-(u.^2)/2) .* (u.^(D-1)-1) * (2*pi).^(-(D+1)/2) .* S ./ (W.^D);
  %
  % KJW's stat_volume can be run like  
  % [PEAK_THRESHOLD, EXTENT_THRESHOLD, PEAK_THRESHOLD_1 EXTENT_THRESHOLD_1] ...
  %  = stat_threshold(S,S/VoxVolMM3,fwhm,Inf,DoesNotMatter,u,.05);
  % EXTENT_THRESHOLD is the critical cluster size in mm3 for clusterp<.05
  %
  % In FSL implmentation, the smoothness file has a DLH value which
  % is in voxel-style units. The fwhm (in mm) can be computed with
  % fwhm_mm = VoxSizeMM*sqrt(4*log(2))/(dLh.^(1/3));
  % dLh = ((fwhm_mm/VoxSizeMM)/sqrt(4*log(2))).^-3
  %Em = exp(-(u.^2)/2) .* (u.^(D-1)-1) * (2*pi).^(-(D+1)/2) .* S .* dLh;
end

% Equation 3
beta = (gamma(D/2+1).*Em./(S.*phiu)).^(2/D);

% Prob that number of voxels in a cluster (n) exceeds k (Bet Eq 2 and 3)
Pnk = exp(-beta.*(k.^(2/D)));
%fprintf('Pnk = %g\n',Pnk);

% Prob of cluster of size k or greater
cdf = 1 - exp(-Em.*Pnk);

return;



