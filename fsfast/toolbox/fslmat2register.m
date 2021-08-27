function R = fslmat2register(Mfsl,cres,rres,sres,nc,nr,ns)
% R = fslmat2register(Mfsl,cres,rres,sres,nc,nr,ns)
%
% Mfsl is the matrix created by FSL-FLIRT with a COR as
% the reference. The object volume has voxel size
% cres X rres X sres (column, row, and slice resolution
% in mm) and is of dimension nr X nc X ns.


%
% fslmat2register.m
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

R = [];

if(nargin ~= 7)
  fprintf('R = fslmat2register(Mfsl,cres,rres,sres,nc,nr,ns)\n');
  return;
end

% FSL D matrix for the COR (voxel size always 1x1x1).
Dc = eye(4);

% FSL D matrix for functional volume %
Df = diag([cres rres sres 1]);

% MGH-FreeSurfer Matrix for COR
Fc = [-1  0  0  (256-1)/2;
       0  0  1 -(256-1)/2;
       0 -1  0  (256-1)/2;
       0  0  0      1];

% MGH-FreeSurfer 
Ff = [-cres  0    0   cres*(nc-1)/2
       0     0  sres -sres*(ns-1)/2
       0   -rres  0   rres*(nr-1)/2
       0     0    0      1];


R = Ff*inv(Df)*inv(Mfsl)*Dc*inv(Fc);

return;
