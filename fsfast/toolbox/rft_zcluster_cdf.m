function cdf = rft_zcluster_cdf(csize,zthresh,fwhm,ssize,D)
% cdf = rft_zcluster_cdf(csize,zthresh,fwhm,ssize,D)
%
% Prob that a cluster >= csize defined by threshold zthresh will be
% found in a D-dim z-field with fwhm smoothness of ssize search space.
% zthresh is the voxel-wise z thresh. ssize and csize are measured
% in non-resel units. csize, fwhm, and ssize are measured in the
% same units.
%
% Based on:
% Friston, Worsley, Frackowiak, Mazziotta, Evans. Assessing the
% significance of focal activations using their spatial extent. HBM
% 1994, 1:214-220. 
%
% Also in:
% Based on Friston, Holmes, Poline, Price, and Frith. Detecting
% Activations in PET and fMRI: Levels of Inference and Power.
% Neuroimage 40, 223-235 (1996).
% 
% $Id: rft_zcluster_cdf.m,v 1.6 2011/03/02 00:04:07 nicks Exp $

%
% rft_zcluster_cdf.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:07 $
%    $Revision: 1.6 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
%if(D==3)
% The only difference here is (u.^(D-1)-1), which appears to come
% from Worsely 1996.
%  Em = exp(-(u.^2)/2) .* (u.^(D-1)-1) * (2*pi).^(-(D+1)/2) .* S ./ (W.^D);
%end

%rhoDng = exp(-(u.^2)/2) .* (u.^(D-1) - 1) * (2*pi).^(-(D+1)/2);
%fprintf('Em = %g\n',Em);
%fprintf('rhoDng = %g\n',rhoDng);
%keyboard

% Equation 3
beta = (gamma(D/2+1).*Em./(S.*phiu)).^(2/D);

% Prob that number of voxels in a cluster (n) exceeds k (Bet Eq 2 and 3)
Pnk = exp(-beta.*(k.^(2/D)));
%fprintf('Pnk = %g\n',Pnk);

% Prob of cluster of size k
cdf = 1 - exp(-Em.*Pnk);

return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These were old attempts. Above seems to disagree with the FSL
% implementation but agrees with KJW's stat_threshold.m and seems
% to agree with actual simulations.

W = fwhm/sqrt(4*log(2));
dLh = W.^(-D); % Same as FSL?
%dLh = 0.493971;

if(1)
Em = ssize .* ((2*pi).^(-(D+1)/2)) .* dLh .* (zthresh.^(D-1)) .* ...
     exp(-(zthresh.^2)/2); 
else
% This appears to be the FSL equation from the techrep
Em = ssize .* ((2*pi).^(-(D+1)/2)) .* dLh .* (zthresh.^(D-1) - 1) .* ...
     exp(-(zthresh.^2)/2); 
end

beta = ( ((gamma(D/2+1) .* Em)) ./ (ssize.*mri_zcdf(zthresh))).^(2./D);

% Prob than n >= k, ie, the number of voxels in a cluster >= csize
Pnk = exp(-beta.*(csize.^(2./D)));
cdf = 1 - exp(-Em.*Pnk);

%------------------------------------------------------------------
%zthresh = fast_p2z(zthresh);
zthresh = zthresh;
pthresh = fast_z2p(zthresh);
rhoD = zthresh .* exp(-(zthresh.^2)/2) / ((2*pi)^1.5); % rho2 for surf
resels = ssize./(fwhm.^D);
invol = resels .* (4*log(2)).^(D/2);
EL = invol.*rhoD;
cons = (((gamma(D/2+1)*((4*log(2))^(D/2))/(fwhm^D)))*rhoD)/pthresh;
pS = exp(-(csize*cons).^(2/D));
P_val_extent = 1-exp(-pS*EL);
cdf2 = P_val_extent;

keyboard

%------------------------------------------------------------------

%keyboard

return;


