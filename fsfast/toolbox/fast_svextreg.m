function err = fast_svextreg(extreg,stem)
% err = fast_svextreg(extreg,stem)
%
% Saves a matrix in a way that it can be used
% as an external regressor for FS-FAST.
%
% extreg is the nframes-by-nregressors matrix
% stem is the bfloat stem to save it in
%
% Eg,
%   nframes = 100; % timepoints
%   nreg    =   3; % number of regressors
%   extreg  = randn(nframes,nreg);
%   fast_svextreg(extreg,'mystem');
%
% This will create mystem_000.bfloat, mystem_000.hdr, and
% mystem.bhdr. When running mkanalysis-sess, use -extreg mystem.
%   
%


%
% fast_svextreg.m
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

if(nargin ~= 2)
  fprintf('err = fast_svextreg(extreg,stem)\n');
  return;
end

[nframes nreg ] = size(extreg);

bmri.te = 0;
bmri.tr = 0;
bmri.ti = 0;
bmri.flip_angle = 0;
bmri.voldim = [nreg 1 1];
bmri.nframes = nframes;
bmri.T = eye(4);
bmri.volres = [1 1 1];
d = fast_mat2vol(extreg,[nreg 1 1]);
err = fast_svbslice(d,stem,[],'bfloat',bmri);

return;
