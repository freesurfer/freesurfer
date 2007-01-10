function sacf = fast_spatialacf(fslice,r0,c0,maxlag,normflag)
% sacf = fast_spatialacf(fslice,r0,c0,maxlag,<normflag>)
% 
% Computes the spatial (2D) autocorrelation function. The
% correlation coefficient is computed as the normalized cross
% correlation of the time courses.
%
% fslice - functional slice [nr nc nf] = size(fslice)
%
% If normflag then normalizes. Otherwise assumes that the time
% courses are already normalized.
%
% gauss std = fwhm/2.36
%
%
%
% (c) Douglas N. Greve, 2004.


%
% fast_spatialacf.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:32 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%

if(nargin < 4 | nargin > 5)
  fprintf('sacf = fast_spatialacf(fslice,r0,c0,maxlag,<normflag>)\n');
  return;
end

[nr nc nf] = size(fslice);
if(r0-maxlag < 1 | r0+maxlag > nr)
  fprintf('ERROR: r0/max lag\n');
  return;
end
if(c0-maxlag < 1 | c0+maxlag > nc)
  fprintf('ERROR: r0/max lag\n');
  return;
end

% 2D sptial acf
sacf = zeros(2*maxlag+1);

v0 = squeeze(fslice(r0,c0,:));
if(normflag) v0 = v0/sqrt(sum(v0.^2)); end

nthdr = 1;
for dr = -maxlag : maxlag
  r = r0 + dr;
  nthdc = 1;
  for dc = -maxlag : maxlag
    c = c0 + dc;
    v = squeeze(fslice(r,c,:));
    if(normflag) v = v/sqrt(sum(v.^2)); end
    sacf(nthdr,nthdc) = v'*v0;
    nthdc = nthdc + 1;
  end
  nthdr = nthdr + 1;
end

return;
