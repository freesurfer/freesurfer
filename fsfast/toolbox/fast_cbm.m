function [cbm, idealcbm, cbmerr] = fast_cbm(seq,maxorder)
% [cbm idealcbm cbmerr] = fast_cbm(seq,maxorder)
%
% Defines nth order counter balancing as the probability of 
% getting a particular stimulus n stimuli after another 
% stimulus. The counter-balancing matrix (CBM) is NxNxmaxorder
% (N = number of event types) where the entry in the ith row 
% and jth column of the kth matrix is the probability of seeing 
% event j k events after event i. 
%
% All the numbers within a column of the idealcbm are the same 
% and equal the number of times the jth type was presented divided 
% by the total number of presentations. The idealcbm is NxN and
% is the same for all orders.
%
% The cbmerr is the difference between the actual cbm and the 
% idealcbm. 
% 
% Notes: event type 0 is removed before processing.
%
% See also fast_cbm2(). Unlike fast_cbm2(), the CBM of all orders
% (upto and including maxorder) are computed.
%


%
% fast_cbm.m
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

cbm = [];
idealcbm = [];
cberr = [];
cbmerr = [];

if(nargin ~= 2)
  fprintf('USAGE: fast_cbm(seq,maxorder)');
  return;
end

% remove 0s from seq %
ind = find(seq ~= 0);
seq = seq(ind);

idlist = unique(seq);
nids = length(idlist);

% Total number of presentations 
nprestot = length(seq);

% Number of times each id was presentated
npresid = zeros(1,nids);
for id = 1:nids
  npresid(id) = length(find(seq==id));
end

% Probability that an id was presented %
ppresid = npresid/nprestot;

idealcbm = repmat(ppresid, [nids 1]);

ncbm = zeros(nids,nids,maxorder);

for order = 1:maxorder,
  for n = 1:(nprestot-order),
    id1 = seq(n);
    id2 = seq(n+order);
    ncbm(id1,id2,order) = ncbm(id1,id2,order) + 1;
  end
end

cbm = ncbm ./ repmat(npresid,[nids 1 maxorder]);

cbmerr = (cbm - repmat(idealcbm,[1 1 maxorder]))./repmat(idealcbm,[1 1 maxorder]);

return;
