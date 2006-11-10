function [hbar, htxt] = barsegstats(segnames,segstats)
%  [hbar htxt] = barsegstats(segnames,segstats)
%
% [segnames segindex segstats] = load_segstats('aseg.stats','bert');
% [hbar htxt] = barsegstats(segnames,segstats(:,3))
% 
% 
% $Id: barsegstats.m,v 1.1 2006/11/10 07:31:09 greve Exp $

mn = min(segstats(:));
mx = max(segstats(:));
d = mx-mn;
a = mn+.02*d;

nseg = size(segnames,1);
x = 1:nseg;
hbar = bar(segstats,'r');
for n = 1:nseg
  htxt(n) = text(n+.25,a,segnames(n,:));
  set(htxt(n),'rotation',90);
  set(htxt(n),'fontsize',4)
end


return;
