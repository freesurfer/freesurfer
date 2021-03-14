function [ovlkeep, ikeep, ovlcmap, cscale] = fast_mkoverlay(ovl,ovlthresh,ovlsat,ovltail,ovlcmap)
% [ovlkeep, ikeep, ovlcmap, cscale] = fast_mkoverlay(ovl,ovlthresh,ovlsat,ovltail,ovlcmap)
%
% See fast_overlay for more info.


%
% fast_mkoverlay.m
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


ovlkeep = [];
ikeep = [];
cscale = [];

if(nargin < 3 | nargin > 5)
  fprintf('USAGE: [ovlkeep, ikeep, ovlcmap, cscale] = fast_mkoverlay(ovl,ovlthresh,ovlsat,<ovltail>,<ovlcmap>)\n');
     return;
end

% Apply defaults %
if(~exist('ovltail')) ovltail = 'posneg'; end
if(isempty(ovltail))  ovltail = 'posneg'; end
if(~exist('ovlcmap')) ovlcmap = nmrcolormap(64,'posneg'); end
if(isempty(ovlcmap)) ovlcmap = nmrcolormap(64,'posneg'); end

% Get lengths of the color map, make sure overlay color map has an 
% even number of entries
novlcmap = size(ovlcmap,1);
if(mod(novlcmap,2) ~= 0)
  fprintf('ERROR: length of ovl cmap must be even\n');
  return;
end

% Compute the scale for the overlay color map %
ovlrange = ovlsat-ovlthresh;
r = [0:novlcmap/2-1]';%'
r = r/max(r);
cscale0 = ovlthresh + ovlrange*r; 
cscale = [-flipud(cscale0); cscale0] ;

% Reshape to 1D 
ovl  = reshape1d(ovl);

% Extract voxels that meet minimum threshold %
switch (ovltail)
 case 'pos', 
   ikeep = find(ovl > ovlthresh);
 case 'neg', 
   ikeep = find(ovl < -ovlthresh);
 case 'posneg', 
   ikeep = find(abs(ovl) > ovlthresh);
 case 'abs', 
   ovl = abs(ovl);
   ikeep = find(ovl > ovlthresh);
 otherwise,
   fprintf('ERROR: tail type %s unrecognized\n',ovltail);
   return;
end

% Make sure there are voxels that meet threshold. If not,
% just return.
nkeep = length(ikeep);
if(nkeep == 0) return; end

% Only keep the ones that meet thresh for further processing. %
ovlkeep = ovl(ikeep);

% Apply saturation %
switch (ovltail)
 case {'pos','abs'}
   isat  = find(ovlkeep > ovlsat);
   ovlkeep(isat) = ovlsat;
 case 'neg', 
   isat  = find(ovlkeep < -ovlsat);
   ovlkeep(isat) = -ovlsat;
 case 'posneg', 
   isat  = find(ovlkeep > ovlsat);
   ovlkeep(isat) = ovlsat;
   isat  = find(ovlkeep < -ovlsat);
   ovlkeep(isat) = -ovlsat;
 otherwise,
   fprintf('ERROR: tail type %s unrecognized\n',ovltail);
   return;
end

ovlkeep0 = ovlkeep; % Backup

% Determine which are positive and negative
ipos = find(ovlkeep > 0);
ineg = find(ovlkeep < 0);

% Scale values so that they are between 1 and novlcmap/2, 
% making sure to handle the sign appropriately. 
newrange = novlcmap/2 - 1;
if(~isempty(ipos))
 ovlkeep(ipos) = round((ovlkeep0(ipos)-ovlthresh)*(newrange/ovlrange) + 1);
 ovlkeep(ipos) = novlcmap/2 + ovlkeep(ipos);
end
if(~isempty(ineg))
  ovlkeep(ineg) = round((-ovlkeep0(ineg)-ovlthresh)*(newrange/ovlrange) + 1);
  ovlkeep(ineg) = novlcmap/2 + 1 - ovlkeep(ineg);
end

return;


