function r = fast_maskvol(volid,maskid,maskedvol)


%
% fast_maskvol.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:31 $
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

r = 1;

if(nargin ~= 3)
  fprintf('USAGE: fast_maskvol(volid,maskid,maskedvol)\n');
  return;
end

fprintf('Loading %s\n',volid);
vol = fmri_ldbvolume(volid);
if(isempty(vol))
  fprintf('ERROR: could not load %s\n',volid);
  return;
end

fprintf('Loading %s\n',maskid);
mask = fmri_ldbvolume(maskid);
if(isempty(mask))
  fprintf('ERROR: could not load %s\n',maskid);
  return;
end

szvol = size(vol);
if(length(szvol) == 4) nt = szvol(4);
else                   nt = 1;
end

nv = prod(szvol(1:3));

vol = reshape(vol,[nv nt]);

ind = find(mask < .5);

vol(ind,:) = 0;

vol = reshape(vol,szvol);

fmri_svbvolume(vol,maskedvol);

r = 0;

return;
