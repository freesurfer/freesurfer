function [seq, nswap] = fmri_seqmutate(oldseq,nswap)
%
% [seq nswap] = fmri_seqmutate(oldseq, <nswap>)
%
%
% $Id: fmri_seqmutate.m,v 1.1 2003/03/04 20:47:40 greve Exp $
%
nseq = length(oldseq);
if(nargin == 1) nswap = nseq; end

m = randperm(nseq);
indold = m(1:nswap);

ind = indold(randperm(nswap));

seq = oldseq;
seq(ind) = oldseq(indold);

return;
