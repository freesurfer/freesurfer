function Xfir = fast_st2fir(st,ntp,psdwin,usew)
% Xfir = fast_st2fir(st,ntp,psdwin,usew)
% $Id: fast_st2fir.m,v 1.1 2004/10/16 02:12:54 greve Exp $

Xfir = 0;

if(nargin ~= 4)
  fprintf('Xfir = fast_st2fir(st,ntp,psdwin,usew)\n');
  return;
end

npres = size(st,1);


for nthpres = 1:npres
  tonset   = st(nthpres,1);
  duration = st(nthpres,2);
  weight   = st(nthpres,3);
  psdwinpres = psdwin;
  psdwinpres(3) = psdwin(3) + duration; 
  taxis = fast_psdwin(psdwinpres,'irftaxis') + tonset;
  
  
end




