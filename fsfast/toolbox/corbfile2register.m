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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:29 $
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
