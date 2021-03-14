function R = corbfile2register(cordir,bstem)
% R = corbfile2register(cordir,bstem)
%
% Computes the MGH-style registration matrix given
% the directory to the COR volume and the stem to
% a bshort/bfloat volume. This should give (about)
% the same results as mri_make_register.


%
% corbfile2register.m
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

if(nargin ~= 2)
  fprintf('R = corbfile2register(cordir,bstem)\n');
  return;
end

bhdrfile = sprintf('%s.bhdr',bstem);
[Tb bres bdim] = load_bhdr(bhdrfile);
if(isempty(Tb)) return; end

corinfofile = sprintf('%s/COR-.info',cordir);
Tc = load_corinfo(corinfofile);
if(isempty(Tc)) return; end

% MGH-FreeSurfer Matrix for COR
Fc = [-1  0  0  (256-1)/2;
       0  0  1 -(256-1)/2;
       0 -1  0  (256-1)/2;
       0  0  0      1];

% MGH-FreeSurfer Matrix for bvolume 
cres = bres(1);
rres = bres(2);
sres = bres(3);
nc   = bdim(1);
nr   = bdim(2);
ns   = bdim(3);
Fb = [-cres  0    0   cres*(nc-1)/2
       0     0  sres -sres*(ns-1)/2
       0   -rres  0   rres*(nr-1)/2
       0     0    0      1];


R = Fb*inv(Tb)*Tc*inv(Fc);


return;
