function [F,phi,r] = tdr_uniform_dft1d(nksamples,nrsamples,revflag,phoffset)
% [F phi r] = tdr_uniform_dft1d(nksamples,nrsamples,<revflag>,<phoffset>)
%
% Creates DFT matrix for uniform sampling using 'default' siemens
% phase trajectory that starts at +pi and goes neg. If revflag=1,
% then this is reversed. phoffset is added to the phase.
%
% F will be nksamples by nrsamples and it will convert
%  a vector in recon space to a vector in phase space.
%
% Note that R = inv(F'*F)*F' = F'/nksamples
% 
% See also tdr_uniform_phtraj, tdr_rtraj.
%
%


%
% tdr_uniform_dft1d.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:35 $
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

F = [];
phi = [];
r = [];
if(nargin < 2 | nargin > 4)
  fprintf('F = tdr_uniform_dft1d(nksamples,nrsamples,<revflag>,<phoffset>)\n');
  return;
end

if(~exist('revflag','var')) revflag = []; end
if(isempty(revflag))        revflag = 0; end

if(~exist('phoffset','var')) phoffset = []; end
if(isempty(phoffset))        phoffset = 0; end

% Phase trajectory
phi = tdr_uniform_phtraj(nksamples);
if(revflag) phi = flipud(phi); end
phi = phi + phoffset;

% Reconn 'trajectory'
r = tdr_rtraj(nrsamples);

% Matrix that converts 
%  a vector in recon space to 
%  a vector in phase space
F = exp(-i*phi*r);

return;
