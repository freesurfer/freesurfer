function h = fast_viewvol(vol,volres,frameres)
% h = fast_viewvol(vol,volres,frameres)


%
% fast_viewvol.m
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

h = [];

if(nargin < 1 | nargin > 3)
  fprintf('USAGE: h = fast_viewvol(vol,volres,frameres)\n');
  return;
end

if(nargin < 2)      volres = [1 1 1]; end
if(isempty(volres)) volres = [1 1 1]; end

if(nargin < 3)        frameres = 1; end
if(isempty(frameres)) frameres = 1; end

