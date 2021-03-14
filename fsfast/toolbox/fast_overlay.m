function [baseovl, ovlcmap, cscale] = fast_overlay(base,ovl,ovlthresh,ovlsat,ovltail,ovlcmap,basecmap)
% baseovl = fast_overlay(base,ovl,ovlthresh,ovlsat,<ovltail>,<ovlcmap>,<basecmap>)
%
% Merge a base and overlay into a true color image. To view, run
% image(baseovl). Run colormap(ovlcmap) before trying to view the 
% colorbar. cscale can be used to assign colors to values.
%
% The base and ovl images can be of any dimension, but each image
% must have the same number of elements in each dimension. The resulting
% image will have one more dimension which will have three elements (RGB).
% cscale is the overlay value associated with each each row of the 
% overlay color map.
%
% base is the base image.
% ovl is the image to be thresholded and overlaid on the base
% ovlthresh is the minimum value ovl needs to be seen in the final image
% ovlsat is the saturation value.
% ovltail (optional) indicates the sign to overlay. pos (positive only),
%   neg (negative only), posneg (positive and negative), and abs
%   (take absolute value and treat as pos). If the ovltail argument
%   is not included or is empty, then posneg is assumed.
% ovlcmap - color map to use for the overlay. The number of rows must
%   be even, and the number of columns must be 3. If not present or
%   empty, the nmrcolormap is used. The overlay color map is returned
%   as part of the output. See below for more on the ovleray colormap.
% basecmap - color map to use for the base. The number of columns must be 3. 
%   If not present or empty, the gray color map (64) is used.
%
% The overlay color map is assumed to be divided into a part for 
% positive values and a part for negative values. For the positive 
% part, the positive values will be mapped to colormap indices 
% between novlcmap/2 + 1 to novlcmap, with increasing values mapped 
% to increasing indices. The negative values will be mapped to 
% indices between 1 and novlcmap/2, with 1 being the saturation 
% value and novlcmap/2 being the threshold.
%
%


%
% fast_overlay.m
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

baseovl = [];
cscale = [];

if(nargin < 4 | nargin > 7)
  fprintf('USAGE: baseovl = fast_overlay(base,ovl,ovlthresh,ovlsat,<ovltail>,<ovlcmap>,<basecmap>)\n');
     return;
end

% Make sure dimensions are consistent %
basedim = size(base);
ovldim = size(ovl);
if(length(basedim) ~= length(ovldim))
  fprintf('ERROR: base/overlay dimension mismatch\n');
  return;
end
for n = 1:length(ovldim)
  if(basedim(n) ~= ovldim(n))
    fprintf('ERROR: base/overlay dimension mismatch\n');
    return;
  end
end

% Apply defaults %
if(~exist('ovltail')) ovltail = ''; end
if(~exist('ovlcmap')) ovlcmap = ''; end
if(~exist('basecmap')) basecmap = repmat([0:63]'/63,[1 3]); end %'
if(isempty(basecmap)) basecmap = repmat([0:63]'/63,[1 3]); end %'

% Get lengths of color maps, make sure overlay color map has an 
% even number of entries
nbasecmap = size(basecmap,1);
if(0)%---------------------%
novlcmap = size(ovlcmap,1);
if(mod(novlcmap,2) ~= 0)
  fprintf('ERROR: length of ovl cmap must be even\n');
  return;
end
end %--------------------------------------------------------%

if(0) %--------------------------------------------------------%
% Compute the scale for the overlay color map %
ovlrange = ovlsat-ovlthresh;
r = [0:novlcmap/2-1]';%'
r = r/max(r);
cscale0 = ovlthresh + ovlrange*r; 
cscale = [-flipud(cscale0); cscale0] ;
end %--------------------------------------------------------%

% Reshape everthing to 1D 
base = reshape1d(base);
ovl  = reshape1d(ovl);

% Rescale the base to fit in the base color map
minbase = min(base);
maxbase = max(base);
if(minbase ~= 1 & maxbase ~= nbasecmap)
  baserange = maxbase-minbase;
  base = (base-minbase)*((nbasecmap-1)/baserange) + 1;
end
base = round(reshape1d(base));

% Compute the True Color image of the base 
baseovl = basecmap(base,:);

% Get the overlay
[ovlkeep ikeep ovlcmap cscale] = ...
  fast_mkoverlay(ovl,ovlthresh,ovlsat,ovltail,ovlcmap);

baseovl(ikeep,:) = ovlcmap(ovlkeep,:);
baseovl = reshape(baseovl,[basedim 3]);

return;
