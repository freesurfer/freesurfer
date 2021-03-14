function dot = fmri_lddot(dotfile)
% dot = fmri_lddot(dotfile)


%
% fmri_lddot.m
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

if(nargin ~= 1)
  msg = 'USAGE: dot = fmri_lddot(dotfile)';
  qoe(msg);error(msg);
end

dot = load(dotfile);
[n1 n2] = size(dot);
ntp = n1/2;
nsensors = n2;

dot = reshape(dot, [ntp 2 nsensors]);
dot = permute(dot, [2 3 1]);

return;

