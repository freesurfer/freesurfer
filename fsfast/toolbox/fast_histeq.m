function [edge, bincenter, binmap] = fast_histeq(y,nbins)
% [edge, bincenter, binmap] = fast_histeq(y,nbins)
%
% computes the bin edges that will result in an equal number of
% samples of y in each bin. 
%
% The size of edge will be nbins+1, where the first bin will be
% between edge(1) and edge(2), etc.
%
%
% To check:
% y = randn(10000,1);
% edge = fast_histeq(y,100);
% nk = histc(y,edge);
% plot(nk(1:end-1));
% plot should have approx 1000 = 10000/nbins at each entry
%
%


%
% fast_histeq.m
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

edge = [];

if(nargin ~= 2 & nargin ~= 3)
  fprintf('[edge bincenter binmap] = fast_histeq(y,nbins)\n');
  return;
end

y = y(:);
ysorted = sort(y);
ny = length(y);
% This is a little bit of a hack for cases where there
% cannot be an equal number per bin.
nperbin = round(ny/nbins);                                          
indedge = [1:nperbin:ny];
n0 = length(indedge);                                                   
if(n0 ~= nbins+1)                                                   
  indedge = [indedge ny];                                                   
end                                                                 
if(indedge(end) ~= ny)                                                  
  indedge(end) = ny;                                                    
end                                                                 
edge = ysorted(indedge)';
bincenter = (edge(1:end-1)+edge(2:end))/2;
binmap = zeros(ny,1);
for nthbin = 1:nbins
  if(nthbin == 1)
    ind = find(y < edge(nthbin+1));
  elseif(nthbin == nbins)
    ind = find(y >= edge(nthbin));
  else
    ind = find(y >= edge(nthbin) & y < edge(nthbin+1));
  end
  binmap(ind) = nthbin;
end

return;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This stuff was not robust 

if(exist('r') ~= 1) r = []; end
if(isempty(r)) r = 20; end

ny = length(y);
nbins0 = r*nbins;

if(ny < nbins0)
  fprintf('ERROR: not enough data points for %d bins (r=%d)\n',nbins,r);
  return;
end

[nperbin bincenter] = hist(y,nbins0);
p = nperbin/ny;
cp = cumsum(p);
binwidth = mean(diff(bincenter));

edge = zeros(nbins+1,1);
edge(1)   = bincenter(1) - binwidth/2;
edge(end) = bincenter(end) + binwidth/2;
for n = 1:nbins-1
  pedge = n/nbins;
  [m i] = min(abs(cp-pedge));
  edge(n+1) = bincenter(i) + binwidth/2;
end

if(nargout < 2) return; end
bincenter = (edge(1:end-1)+edge(2:end))/2;

if(nargout < 3) return; end

binmap = zeros(size(y));
for n = 1:nbins
  if(n == 1)
    ind = find(y <= edge(n+1) );
  elseif(n==nbins)
    ind = find(edge(n) < y);
  else
    ind = find(edge(n) < y & y <= edge(n+1) );
  end
  binmap(ind) = n;
end



%nk = histc(y,edge);
%plot(nk(1:end-1));
%keyboard

return;























