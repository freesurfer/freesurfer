function [acfreg, nkeep, pvs, pvsthresh, uacf, usim] = fast_acf_svdreg(acf,pvsthresh)
% [acfreg,nkeep,pvs,u] = fast_acf_svdreg(acf,pvsthresh)


%
% fast_acf_svdreg.m
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

acfreg = [];
uacf = [];
usim = [];

if(nargin ~= 1 & nargin ~= 2)
  fprintf('USAGE: [acfreg,nkeep,pvs,u] = fast_acf_svdreg(acf,pvsthresh)\n');
  return;
end

[nf nv] = size(acf);

if(nargin == 1)
  % Do simulation with white noise %
  [pvsthresh usim] = fast_sim_acfsvd(nf, 2000);
else
  if(nf-1 ~= length(pvsthresh))
    fprintf('ERROR: nf-1 ~= length(pvsthresh)\n');
    return;
  end
end

acf = acf(2:nf,:);
Macf = acf*acf'; %'
[uacf s v] = svd(Macf);
ds = diag(s);
pvs = 100*ds/sum(ds);

%figure(2);
%plot(1:20,pvs(1:20),1:20,pvsw(1:20));
%keyboard

for nkeep = 1:nf-1
  if(pvs(nkeep) < pvsthresh(nkeep)) break;end
end
nkeep = nkeep - 1;

if(nkeep ==0)
  acfreg = zeros(nf,nv);
  acfreg(1,:) = 1;
end

% nkeep = max(find(pvs > pvsw));
uwht = uacf(:,nkeep+1:nf-1);
acfreg = (eye(nf-1)-uwht*uwht')*acf; %'
acfreg = [ones(1,nv); acfreg];

return;
