function Cnorm = fast_norm_con(C)
% Cnorm = fast_norm_con(C)
%
% Make sure that the positives of each row sum to 1 and that the
% negatives of each row sum to -1. This also assures that each row
% sums to zero if there are positives and negatives in the row.
%
% This is a replicated of code found in fast_contrastmtx.m
%

Cnorm = [];
if(nargin ~= 1)
  fprintf('Cnorm = fast_norm_con(C)\n');
  return;
end

Cnorm = C;
for nthrow = 1:size(C,1);
  % positives %
  ind = find(C(nthrow,:)>0);
  if(~isempty(ind))
    xsum = sum(C(nthrow,ind));
    Cnorm(nthrow,ind) = C(nthrow,ind)/xsum;
  end
  % negatives %
  ind = find(C(nthrow,:)<0);
  if(~isempty(ind))
    xsum = sum(C(nthrow,ind));
    Cnorm(nthrow,ind) = C(nthrow,ind)/abs(xsum);
  end
end

return;
