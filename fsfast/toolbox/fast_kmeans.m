function [kmeans, kmap, dmin, niters, yhat] = fast_kmeans(y,nc,kmeans0,nitersmax,nfix,ndiffmax)
% [kmeans, kmap, dmin, niters, yhat] = fast_kmeans(y,nc,<kmeans0>,<nitersmax>,<nfix>,<ndiffmax>)
%
% y is n-variates by nsamples
% nc is number of classes (a better name would have been nk)
% If nitersmax is not specified, uses 100.
% If kmeans0 is not specified, uses first nc of y:
%      kmeans0 = y(:,1:nc);
% The mean error is mean(dmin) = mean(abs(y-yhat)). Note
% that this is an L1, not L2, measure.
%
% nfix - fix the first nfix class means as specified in
% kmeans0. nc-nfix class means are adapted.
%
% kmeans is n-variates by n-classes
% kmap is 1 by nsamples
%

%
% fast_kmeans.m
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

kmeans = [];
kmap = [];
dmin = [];
niters = 0;

if(nargin < 2 | nargin > 6)
 fprintf('[kmeans kmap dmin niters yhat] = fast_kmeans(y,nc,<kmeans0>,<nitersmax>,<nfix>,<ndiffmax>)\n');
 return;
end

if(~exist('nfix','var')) nfix = []; end
if(isempty(nfix)) nfix = 0; end
if(~exist('kmeans0','var') & nfix > 0)
  fprintf('ERROR: must specify kmeans0 with nfix\n');
  return;
end
if(nfix >= nc)
  fprintf('ERROR: nfix >= number of classesx\n');
  return;
end

if(~exist('ndiffmax','var')) ndiffmax = []; end
if(isempty(ndiffmax)) ndiffmax = 0; end

if(~exist('nitersmax','var')) nitersmax = []; end
if(isempty(nitersmax))        nitersmax = 100; end

if(~exist('kmeans0','var'))   kmeans0 = []; end

% If unspecified, init kmeans to first nc examples
if(isempty(kmeans0))  kmeans0 = y(:,1:nc); end

if(nc ~= size(kmeans0,2))
  fprintf('ERROR: nc does not equal kmean0\n');
  return;
end

[nf nv] = size(y);
%tic;
niters = 0;
ndiff = nv;
kmeans = kmeans0;
while(niters < nitersmax & ndiff > ndiffmax)

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
    if(~isempty(ind)) 
      kmeans(:,c) = mean(y(:,ind),2);
    else
      % This class did not come closest to any samples
      % Assign it to the sample that it is closest to
      % This used to assign to some samples randomly
      % which created non-deterministicness
      [mm ii] = min(d(c,:));
      kmeans(:,c) = y(:,ii);
      kmap(ii) = c;
    end
  end

  % Check whether anything has changed 
  if(niters ~= 0 )
    ndiff = length(find( (kmap-kmap0) ~= 0));
  end

  kmap0 = kmap;
  niters = niters + 1;
  if(mod(niters,10)==0 | niters == 1)
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




