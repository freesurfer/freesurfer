function [havg, errstd] = fmri_untangle2(havgstd, hdrdat)
% 
% [havg, errstd] = fmri_untangle2(havgstd, hdrdat)
%
if(nargin ~= 2)
  msg = 'Usage: [havg errstd] = fmri_untangle2(havgstd, hdrdat)';
  qoe(msg);error(msg);
end

Nh   = hdrdat.Nh;
Nc   = hdrdat.Nc;
Nnnc = hdrdat.Nnnc;

tmp1 = repmat([1:Nh]', [1 Nnnc]); %'
tmp2 = repmat([2*Nh : 2*Nh : 2*Nh*Nc-1], [Nh 1]);
indavg = reshape1d(tmp1+tmp2);

sz = size(havgstd);
dim = length(sz);
if(dim==1)
  havg = havgstd(indavg);
  errstd = havgstd(Nh+1);
  return;
end

szv = sz(1:dim-1);
nv = prod(szv);
nhtot = sz(dim);

%keyboard
havgstd = reshape(havgstd, [nv nhtot]);

havg = havgstd(:,indavg);
errstd = havgstd(:,Nh+1);


havg   = reshape(havg, [szv size(havg,2)]);
errstd = reshape(errstd, [szv 1]);

return;

