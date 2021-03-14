function [cbm, idealcbm] = fast_cbm2(seq,order)
%
% [cbm idealcbm] = fast_cbm2(seq,order)
%
% Defines nth order counter balancing as the probability of 
% getting a particular stimulus after a order-length sequence
% of stimuli. If there are N event types in the sequence, then
% the number of possible order-length sequnces is M = N^order
% The counter-balance matrix (CBM) has a size MxN. The entry in 
% the ith row and jth column is the probability of seeing event
% type j after subsequence i.  
%
% Eg. N=3, order=3. M = 27. CBM is MxN:
% i=1 --> (1,1,1)
% i=2 --> (1,1,2)
% i=3 --> (1,1,3)
% i=4 --> (1,2,1)
% i=5 --> (1,2,2)
% ...
%
% All the numbers within a column of the idealcbm are the same 
% and equal the number of times the jth type was presented divided 
% by the total number of presentations. The idealcbm is also MxN.
%
% See also fast_cbm(). Unlike fast_cbm(), the CBM of only the 
% specified order is computed. 
%
% Order must be <= 3.


%
% fast_cbm2.m
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

if(nargin ~= 2)
  fprintf('USAGE: fast_cbm2(seq,order)');
  return;
end

if(order < 1 | order > 3)
  fprintf('ERROR: order = %d, range: 1 <= order <=3 \n',order);
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

sz = nids*ones(1,order);
switch(order)
  case 1, 
    midlist = sub2ind(sz,seq(1:nprestot-1));
  case 2, 
    midlist = sub2ind(sz,seq(1:nprestot-2),seq(2:nprestot-1));
  case 3, 
    midlist = sub2ind(sz,seq(1:nprestot-3),seq(2:nprestot-2),seq(3:nprestot-1));
end

nmids = nids.^order;
ncbm = zeros(nmids,nids);
for n = order+1:nprestot,
  mid = midlist(n-order);
  id = seq(n);
  ncbm(mid,id) = ncbm(mid,id) + 1;
end

idealcbm = repmat(ppresid, [nmids 1]);

cbm = ncbm ./ repmat(npresid,[nmids 1]);

return;
