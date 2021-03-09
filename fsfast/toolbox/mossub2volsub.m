function [rv, cv, sv, im, tszmos] = mossub2volsub(rm, cm, szvol, tszmos)
% [rv cv sv tszmos] = mossub2volsub(rm, cm, szvol, tszmos)
%
% Computes the subscripts in a volume (row, col, slice) that correspond
% to a subscript in a mosaic (row, col).  Note that when the mosaic has
% been padded with blank images, the subscripts in those images do not
% have corresponding subscripts in the volume.
%
% rm - row in the mosaic
% cm - column in the mosaic
% szvol - size of the volume (Nrows, Ncols, Nslices, ...)
% tszmos - size (rows, cols) of the mosaic measured in tiles (optional)
%
% rv - row in the volume
% cv - column in the volume
% sv - slice in the volume
%
% If tszmos is not specified, a default one will be computed using
% the function defmossize.
%
% See also: mos2vol vol2mos mosind2volind mossub2volsub 
%           volind2mosind volsub2mossub defmossize


%
% mossub2volsub.m
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

if(nargin ~= 3 & nargin ~= 4)
  msg = 'USAGE: [rv cv sv tszmos] = mossub2volsub(rm, cm, szvol, <tszmos>)';
  error(msg);
end

if(length(rm) ~= length(cm))
  msg = sprintf('rm (%d) and cm (%d) do not have the same length',...
               length(rm),length(cm));
  error(msg);
end

Nvr = szvol(1);
Nvc = szvol(2);
Nvs = szvol(3);

if(nargin == 3) tszmos = []; end
tszmos = defmossize(Nvs, tszmos);

Ntr = tszmos(1);
Ntc = tszmos(2);

rt = floor((rm-1)/Nvr) + 1; % tile row in mosaic
ct = floor((cm-1)/Nvc) + 1; % tile col in mosaic

%fprintf('rt = %3d, ct = %3d\n',rt,ct);

rv = rm - (rt-1)*Nvr;
cv = cm - (ct-1)*Nvc;
sv = ct + (rt-1)*Ntc;

% Exclude out of range slices. Rows and cols should always
% be in range, but slices can be out of range when there are more
% mosaic tiles than volume slices.
im = find(sv <= Nvs);
if(length(im) ~= length(sv))
  rv = rv(im);
  cv = cv(im);
  sv = sv(im);
end

return
