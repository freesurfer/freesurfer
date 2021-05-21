function [seqkeep, trkeep, travg, trstd] = fmri_minseqtr(seqseed,Nh,nsearch,nkeep,nswap,Mss)
%
% [seqkeep, trkeep, travg, trstd] = fmri_minseqtr(seqseed,Nh,nsearch,nkeep,nswap,<Mss>)
%
%
%


%
% fmri_minseqtr.m
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

if(nargin ~= 4 & nargin ~= 5 & nargin ~= 6)
  msg = 'USAGE: fmri_minseqtr(seqseed,Nh,nsearch,nkeep,<<nswap>,Mss>)';
  qoe(msg);error(msg);
end

seq = seqseed;
nseq = length(seq);
trmin = 10^10;
trsum = 0;
trsumsq = 0;

seqkeep = zeros(nseq,nkeep);
trkeep  = 10^10*ones(nkeep,1);
if(nargin >= 6) 
  m2 = Mss'*Mss;
  indm = find(diag(m2)==1);
end


for n = 1:nsearch,

  xtx = 0;
  kmax = 10;
  k = 0;
  while(cond(xtx)>1000 & k < kmax)
    if(nargin == 4) seq = fmri_seqmutate(seq);
    else            seq = fmri_seqmutate(seq,nswap);
    end

    x = fmri_seq2scm(seq,Nh);

    if(nargin < 6)     xtx = x'*x;
    else               xtx = x(indm,:)' * x(indm,:);
    end
    %fprintf('%2d %g\n',k,cond(xtx));
    k = k + 1;
  end
  if(k >= kmax)
    msg = sprintf('Cound not find xtx with condition < 1000');
    qoe(msg);error(msg);
  end

  tr = trace(inv(xtx));

  if(tr < trkeep(nkeep))
    trkeep(nkeep) = tr;
    seqkeep(:,nkeep) = seq;

    [tmp1 indx] = sort(trkeep);
    trkeep = tmp1;
    seqkeep = seqkeep(:,indx);
  end

  trsum   = trsum + tr;
  trsumsq = trsumsq + tr*tr;

end

travg = trsum/nsearch;
trstd = sqrt(trsumsq/(nsearch-1) - travg*travg);

return;
