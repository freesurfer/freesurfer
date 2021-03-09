function vol = MRIvol2vol(mov,targ,R)
% vol = MRIvol2vol(mov,targ,<R>)
%
% mov  = MRIread('mov'); % only one frame
% targ = MRIread('targ');
%
% R is tkregister2 matrix (maps targ-to-mov)
%   If unspec, uses header reg based on 
%   scanner vox2ras.
%
% Currently only uses nearest neighbor.
%
%


%
% MRIvol2vol.m
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

vol = [];
if(nargin < 2 | nargin > 3)
  fprintf('vol = MRIvol2vol(mov,targ,<R>)\n');
  return;
end

Sm = mov.vox2ras0;
Tm = mov.tkrvox2ras;
St = targ.vox2ras0;
Tt = targ.tkrvox2ras;

if(nargin == 2)  R = Tm*inv(Sm)*St*inv(Tt); end
  
% Target vox to Mov vox Matrix
Vt2m = inv(Tm)*R*Tt;

nct = targ.volsize(2);
nrt = targ.volsize(1);
nst = targ.volsize(3);
nvt = prod(targ.volsize);
[tc tr ts] = meshgrid([0:nct-1],[0:nrt-1],[0:nst-1]);
tcrs = [tc(:) tr(:) ts(:) ones(nvt,1)]';

fprintf('Computing indices ... ');tic;
mcrs = round(Vt2m * tcrs);
fprintf(' ... done %g\n',toc);

ncm = mov.volsize(2);
nrm = mov.volsize(1);
nsm = mov.volsize(3);
nvm = prod(mov.volsize);

fprintf('Getting ok ... ');tic
mc = mcrs(1,:);
mr = mcrs(2,:);
ms = mcrs(3,:);
indok = find(mc >= 0 & mc < ncm & ...
	     mr >= 0 & mr < nrm & ...
	     ms >= 0 & ms < nsm);
fprintf(' ... done %g\n',toc);
nok = length(indok);
fprintf('nok = %d\n',nok);

fprintf('Getting tind ... ');tic
tc = tc(indok);
tr = tr(indok);
ts = ts(indok);
tind = sub2ind(targ.volsize,tr+1,tc+1,ts+1);
fprintf(' ... done %g\n',toc);

fprintf('Getting mind ... ');tic
mc = mc(indok);
mr = mr(indok);
ms = ms(indok);
mind = sub2ind(mov.volsize,mr+1,mc+1,ms+1);
fprintf(' ... done %g\n',toc);

fprintf('Resampling ... ');tic
vol = targ;
vol.vol = zeros(nrt,nct,nst,1);
vol.vol(tind) = mov.vol(mind);
fprintf(' ... done %g\n',toc);

return;




