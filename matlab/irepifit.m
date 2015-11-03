function irepistruct = irepifit(irepistruct)
% irepistruct = irepifit(irepistruct)
% Performs LMS fit.
%

s = irepistruct; % copy into s for easy handling

% s = irepitiming(s); Need to do this if fitting seq params

% Need to be able to cache so not recomputing every slice
% Need this way too for incremental
s = irepisynth(s); 

s.X0 = s.yRead;

s.iexclude1 = [];
if(s.nexclude > 0) 
  indSlice = find(s.EventSliceNo == s.sliceno | s.EventSliceNo < 0);
  indROS = find(s.IsReadOut(indSlice));
  nthSliceAcq = s.AcqSliceNo(indSlice(indROS));
  for n = 1:s.nexclude
    kk = find(nthSliceAcq == n);
    s.iexclude1 = [s.iexclude1 kk];
  end
end
s.iexclude2 = [];
if(s.nminexclude > 0)
  [ysort isort] = sort(s.y);
  s.iexclude2 = sort(isort(1:s.nminexclude));
end

if(~isempty(s.iexclude1) & isempty(s.iexclude2))
  s.iexclude = s.iexclude1;
elseif(isempty(s.iexclude1) & ~isempty(s.iexclude2))
  s.iexclude = s.iexclude2;
else
  s.iexclude = unique([s.iexclude1; s.iexclude2]);
end
s.indkeep = setdiff([1:s.ntp],s.iexclude);

s.yFit = s.y(s.indkeep,:);
s.X = s.X0(s.indkeep,:);
s.tFit = s.tRead(s.indkeep);

nX = size(s.X,2);
Xvar = sum(s.X.*s.X); % need to cache, unless exclude
s.M0 = sum(s.X.*repmat(s.yFit,[1 nX]))./Xvar;
s.yhat = s.X.*repmat(s.M0,[length(s.yFit) 1]);
s.res = repmat(s.yFit,[1 nX])-s.yhat;
s.yhat0 = s.X0.*repmat(s.M0,[length(s.y) 1]);
if(0)
  s.M0 = inv(s.X'*s.X)*(s.X'*s.y);
  s.yhat = s.M0*s.X;
  s.res = s.y-s.yhat;
  s.yhat0 = s.M0*s.X0;
end
s.rstd = std(s.res);
err = s.rstd;
s.err = err;

irepistruct = s;

return;
