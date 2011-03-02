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
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:06 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
