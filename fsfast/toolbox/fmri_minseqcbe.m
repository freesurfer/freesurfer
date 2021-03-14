function [seqkeep, cbekeep, cbeavg, cbestd] = fmri_minseqcbe(seqseed,Nh,nsearch,nkeep,idealxtx,nswap)
%
% [seqkeep, cbekeep, cbeavg, cbestd] = 
%    fmri_minseqcbe(seqseed,Nh,nsearch,nkeep,idealxtx,nswap)
%
%
%
%


%
% fmri_minseqcbe.m
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

if(nargin ~= 5 & nargin ~= 6)
  msg = 'USAGE: fmri_minseqcbe(seqseed,Nh,nsearch,nkeep,idealxtx,nswap)';
  qoe(msg);error(msg);
end

seq = seqseed;
nseq = length(seq);
nseq2 = nseq*nseq;

cbesum = 0;
cbesumsq = 0;

seqkeep = zeros(nseq,nkeep);
cbekeep  = 10^10*ones(nkeep,1);

for n = 1:nsearch,

  xtx = 0;
  while(cond(xtx)>1000)
    if(nargin == 5) seq = fmri_seqmutate(seq);
    else            seq = fmri_seqmutate(seq,nswap);
    end
    x = fmri_seq2scm(seq,Nh);
    xtx = x'*x;
  end

  cbe = mean(abs(reshape1d(idealxtx-xtx)));

  if(cbe < cbekeep(nkeep))
    cbekeep(nkeep) = cbe;
    seqkeep(:,nkeep) = seq;

    [tmp1 indx] = sort(cbekeep);
    cbekeep = tmp1;
    seqkeep = seqkeep(:,indx);
  end

  cbesum   = cbesum + cbe;
  cbesumsq = cbesumsq + cbe*cbe;

end

cbeavg = cbesum/nsearch;
cbestd = sqrt(cbesumsq/(nsearch-1) - cbeavg*cbeavg);

return;
