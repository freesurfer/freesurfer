function r = fast_maskvol(volid,maskid,maskedvol)


%
% fast_maskvol.m
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
