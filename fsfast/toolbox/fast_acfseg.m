function [binmap, nperbin] = fast_acfseg(functemplate,nbins,mask)
% [binmap, nperbin] = fast_acfseg(functemplate,nbins,mask)
%
%


%
% fast_acfseg.m
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

ymn = functemplate;
indmask = find(mask);
t0 = toc;
rgpar0  = fast_raygausinit(ymn(indmask));
[rgpar rgcost exitflag] = fminsearch('fast_raygauscost',rgpar0,[],ymn(indmask));
if(exitflag ~= 1)
  fprintf('ERROR: RG optimization, exit=%d\n',exitflag);
  if(~monly) 
    quit; 
  end
  return;
end
fprintf('  done searching (%g min)\n',(toc-t0)/60);
RGalpha = rgpar(1);
Rmu  = rgpar(2);
Gmu  = rgpar(3);
Gstd = rgpar(4);

if(0)
  ythresh = fast_raygausthresh(rgpar);
  [h0 x0] = hist(ymn(indmask),100);
  pdf0 = h0/trapz(x0,h0);
  rgpdf = fast_raygaus(x0,rgpar);
  plot(x0,pdf0,x0,rgpdf);
end

ymax = Gmu + 3*Gstd;
indsupmax = find(ymn > ymax);
nsup = length(indsupmax);
fprintf('RG: alpha=%g, Rmu=%g, Gmu=%g, Gstd=%g, ymax=%g, nsup=%d\n',...
    RGalpha,Rmu,Gmu,Gstd,ymax,nsup);
ymn(indsupmax) = ymax;

% ------- Get the bins of the intensity ----------
fprintf('Computing bins t=%g\n',toc);
[hbin xbin] = hist(ymn(indmask),nbins);
dxbin = xbin(2)-xbin(1);
binmap = zeros(size(ymn));
nperbin = zeros(nbins,1);
for nthbin = 1:nbins
  yminbin = xbin(nthbin)-dxbin/2;
  ymaxbin = xbin(nthbin)+dxbin/2;
  if(nthbin == 1)     yminbin = -inf; end
  if(nthbin == nbins) ymaxbin = +inf; end
  indbin = find(yminbin < ymn(indmask) & ymn(indmask) < ymaxbin);
  nperbin(nthbin) = length(indbin);
  binmap(indmask(indbin)) = nthbin;
  fprintf('bin %2d  %3d  %8.3f %8.3f\n',nthbin,nperbin(nthbin),yminbin,ymaxbin);
end

return;



