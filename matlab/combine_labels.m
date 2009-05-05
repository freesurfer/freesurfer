function lboth = combine_labels(l1, l2)
% lboth = combine_labels(l1, l2)

ind1 = l1(:,1) ;
ind2 = l2(:,1);


same = 0 ;
same1 = zeros(size(ind1)) ;
same2 = zeros(size(ind2)) ;

for i=1:length(ind1)
  for j=1:length(ind2)
    if (ind1(i) == ind2(j))
      same = same+1 ;
      same1(i) = j ;
      same2(j) = i ;
    end
  end
end


for i=1:length(ind1)
  if (same1(i) == 0)
    lboth(i,:) = l1(i,:) ;
  else
    lboth(i,:) = l1(i,:) ;
    lboth(i,5) = l1(i,5) + l2(same1(i),5) ;
  end
end

ind = i ;
for i=1:length(ind2)
  if (same2(i) == 0)
    lboth(ind,:) = l2(i,:) ;
    ind = ind+1 ;
  end
end




