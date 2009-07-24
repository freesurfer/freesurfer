function hs = histo_smooth(h, niters)
% hs = histo_smooth(h, niters)
hsave = h ;
  
if (size(h,1) > 1)
  h = h(:,2); 
end

nbins = length(h);
bins = 1:nbins;
bm1 = [1 1:nbins-1];
bp1 = [2:nbins nbins];

hs = h ;
for n=1:niters
  hs = (hs(bm1) + hs + hs(bp1))/ 3;
end

if (size(hsave) ~= size(h))
  hs1 = hs ;
  hs = hsave ;
  hs(:,2) = hs1 ;
end

