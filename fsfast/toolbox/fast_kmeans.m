function [kmeans, kmap, dmin, niters, yhat] = fast_kmeans(y,nc,kmeans0,nitersmax,nfix)
% [kmeans, kmap, dmin, niters, yhat] = fast_kmeans(y,nc,<kmeans0>,<nitersmax>,<nfix>)
%
% nc is number of classes (a better name would have been nk)
% If nitersmax is not specified, uses 100.
% If kmeans0 is not specified, uses first nc of y.
% The mean error is mean(dmin) = mean(abs(y-yhat)). Note
% that this is an L1, not L2, measure.
%
% nfix - fix the first nfix class means as specified in
% kmeans0. nc-nfix class means are adapted.
%
% $Id: fast_kmeans.m,v 1.6 2004/05/27 01:37:32 greve Exp $
%

kmeans = [];
kmap = [];
dmin = [];
niters = 0;

if(nargin < 2 | nargin > 5)
 fprintf('[kmeans, kmap, dmin, niters, yhat] = fast_kmeans(y,nc,<kmeans0>,<nitersmax>,<nfix>)\n');
 return;
end

if(~exist('nfix','var')) nfix = 0; end
if(isempty(kmeans0) & nfix > 0)
  fprintf('ERROR: must specify kmeans0 with nfix\n');
  return;
end
if(nfix >= nc)
  fprintf('ERROR: nfix >= number of classesx\n');
  return;
end

if(~exist('nitersmax','var')) nitersmax = 100; end
if(~exist('kmeans0','var'))   kmeans0 = []; end

% If unspecified, init kmeans to first nc examples
if(isempty(kmeans0))  kmeans0 = y(:,1:nc); end

if(nc ~= size(kmeans0,2))
  fprintf('ERROR: nc does not equal kmean0\n');
  return;
end

[nf nv] = size(y);

tic;
niters = 0;
ndiff = nv;
kmeans = kmeans0;
while(niters < nitersmax & ndiff ~= 0)

  % Compute the distances between each input and each class
  % Use L1 distance.
  for c=1:nc
    yhatc = repmat(kmeans(:,c),[1 nv]) ;
    r = y - yhatc; % residual
    if(nf > 1)  d(c,:) = mean(abs(r), 1 );
    else        d(c,:) = abs(r);
    end
  end

  % Remap each input to classes
  [dmin kmap] = min(d,[],1);

  % Recompute the class means for unfixed classes
  for c = nfix+1:nc
    ind = find(kmap==c);
    if(~isempty(ind)) kmeans(:,c) = mean(y(:,ind),2);
    else
      % Randomly reassign to one of the cols of y
      ind = round(nv*rand)+1;
      kmeans(:,c) = y(:,ind);
    end
  end

  % Check whether anything has changed 
  if(niters ~= 0 )
    ndiff = length(find( (kmap-kmap0) ~= 0));
  end

  kmap0 = kmap;
  niters = niters + 1;
  if(1 | mod(niters,10)==0 | niters == 1)
    %fprintf('%3d %5d %14.13f %g\n',niters,ndiff,mean(dmin),toc);
  end

end % iteration loop

%fprintf('%3d %5d %14.13f %14.13f %g\n',niters,ndiff,...
%	mean(dmin),sqrt(mean(dmin)),toc);

%------- Create estimate of the input ----------%
if(nargout == 5)
  yhat = zeros(size(y));
  for c = 1:nc
    ind = find(kmap==c);
    nind = length(ind);
    yhat(:,ind) = repmat(kmeans(:,c),[1 nind]);
  end
end

return;
%-------------------------------------------------------------------%
%-------------------------------------------------------------------%
%-------------------------------------------------------------------%

%-------------------------------------------------------------%
% These functions are no longer used, but could come in
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




