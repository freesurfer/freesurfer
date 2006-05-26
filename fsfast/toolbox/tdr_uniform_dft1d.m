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
% $Id: tdr_uniform_dft1d.m,v 1.1 2006/05/26 23:49:11 greve Exp $

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
