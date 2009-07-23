function smooth = histo_smooth(h, niter)
% smooth = histo_smooth(h, niter)

len = length(h) ;
b = [1:len] ;
b1 = [2:len len] ;
b2 = [1 1:len-1] ;

for n=1:niter
  smooth = (h(b) + h(b1) + h(b2)) / 3 ;
  h = smooth ;
end

     

