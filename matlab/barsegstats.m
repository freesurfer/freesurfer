function [hbar, htxt] = barsegstats(segnames,segstats)
%  [hbar htxt] = barsegstats(segnames,segstats)
%
% [segnames segindex segstats] = load_segstats('aseg.stats','bert');
% [hbar htxt] = barsegstats(segnames,segstats(:,3))
% 
% 


%
% barsegstats.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
%    $Revision: 1.4 $
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


mn = min(segstats(:));
mn = 0;
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
