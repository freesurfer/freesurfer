function [kmeans, kmap, d2min, niters] = fast_kmeans(y,nc,kmeans0,nitersmax)
% [kmeans, kmap, d2min, niters] = fast_kmeans(y,nc,kmeans0,nitersmax)

kmeans = [];
kmap = [];
d2min = [];
niters = 0;

if(nargin < 2 | nargin > 4)
 fprintf('[kmeans kmap d2min niters] = fast_kmeans(y,nc,kmeans0,nitersmax)\n');
 return;
end

[nf nv] = size(y);

if(~exist('nitersmax')) nitersmax = 100; end
if(~exist('kmeans0'))   kmeans0 = y(:,1:nc); end

tic;
niters = 0;
ndiff = nv;
kmeans = kmeans0;
while(niters < nitersmax & ndiff ~= 0)

  %d2 = dist2(y,kmeans);
  for c=1:nc
    d2(c,:) = sum( (y - repmat(kmeans(:,c),[1 nv]) ).^2 );
  end

  [d2min kmap] = min(d2,[],1);

  for c=1:nc
    ind = find(kmap==c);
    kmeans(:,c) = mean(y(:,ind),2);
  end

  if(niters ~= 0 )
    ndiff = length(find( (kmap-kmap0) ~= 0));
  end

  kmap0 = kmap;
  niters = niters + 1;
  fprintf('%3d %5d %8.4f %g\n',niters,ndiff,mean(d2min),toc);

end

return;

%-------------------------------------------------------------%
% These functions are nolong used, but could come in
% handy at some point.
%-------------------------------------------------------------%
function d2 = dist2(y,kmeans)
  nv = size(y,2);
  nc = size(kmeans,2);
  for c=1:nc
    d2(c,:) = sum( (y - repmat(kmeans(:,c),[1 nv]) ).^2 );
  end
return
%-------------------------------------------------------------%
function kmeans = compute_kmeans(y,kmap,nc)
  nv = size(y,2);
  for c=1:nc
    ind = find(kmap==c);
    kmeans(:,c) = mean(y(:,ind),2);
  end
return
%-------------------------------------------------------------%




