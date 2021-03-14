function [seq, nswap] = fmri_seqmutate(oldseq,nswap)
%
% [seq nswap] = fmri_seqmutate(oldseq, <nswap>)
%
%
%
%


%
% fmri_seqmutate.m
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


nseq = length(oldseq);
if(nargin == 1) nswap = nseq; end

m = randperm(nseq);
indold = m(1:nswap);

ind = indold(randperm(nswap));

seq = oldseq;
seq(ind) = oldseq(indold);

return;
