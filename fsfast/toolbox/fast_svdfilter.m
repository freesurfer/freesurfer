% fast_svdfilter.m

% instem
% maskstem
% pforder
% xmatfile
% pthresh
% pvsthresh
% 
% outstem
% extregstem

tic

fsd = '/space/greve/1/users/greve/wagner/JL/bold';
funcrstem = 'fmcsm5';
maskrstem = 'brain';
analysis = 'onoff';
dimfrac = .75;
sigthresh = 2;
runidlist = strvcat('006','007','008','009');

pthresh = 10.^(-sigthresh);

maskastem = sprintf('%s/masks/%s',fsd,maskrstem);
mask = fmri_ldbvolume(maskastem);
if(isempty(mask))
  fprintf('ERROR: could not load %s\n',maskastem);
  return; quit;
end
indmask = find(mask);
nmask = length(find(mask));
if(nmask == 0)
  fprintf('ERROR: could not find any voxels in mask\n');
  return; quit;
end
fprintf('INFO: found %d voxels in mask\n',nmask);
  

nruns = size(runidlist,1);
for nthrun = 1:nruns

  runid = deblank(runidlist(nthrun,:));
  fprintf('run %s   (%g) -----------------\n',runid,toc);

  funcastem = sprintf('%s/%s/%s',fsd,runid,funcrstem);
  xmatfile  = sprintf('%s/%s/%s.mat',fsd,runid,analysis);

  XX = load(xmatfile);
  if(isempty(XX))
    fprintf('ERROR: could not load %s\n',xmatfile);
    return; quit;
  end
  X = XX.Xfinal;
  Ntask = XX.Nnnc*XX.Navgs_per_cond;
  Nbeta = size(X,2);
  C = zeros(Ntask,Nbeta);
  C(1:Ntask,1:Ntask) = 1;
  dof = size(X,1)-size(X,2);
  B = inv(X'*X)*X';
  
  fprintf('Loading functionals (%g)\n',toc);
  y = fmri_ldbvolume(funcastem);
  if(isempty(y))
    fprintf('ERROR: could not load %s\n',funcastem);
    return; quit;
  end
  [ns nr nc nf] = size(y);
  nv = ns*nr*nc;
  
  if(nf ~= size(X,1))
    fprintf('ERROR: mismatch in the number of frames\n');
    return; quit;
  end
  if(nv ~= prod(size(mask)))
    fprintf('ERROR: mismatch in the number of voxels\n');
    return; quit;
  end
  y = reshape(y,[nv nf])';
  y = y(:,indmask);
  
  fprintf('Analyzing functionals to get exclusion mask (%g)\n',toc);
  beta = B*y;
  r = y-X*beta;
  rvar = sum(r.^2)/dof;
  [F, dof1, dof2, g] = fast_fratio(beta,X,rvar,C);
  p = FTest(dof1, dof2, F);

  indexclude = find(p < pthresh);
  indinclude = find(p >= pthresh);
  nexclude = length(indexclude);
  ninclude = length(indinclude);
  fprintf('Excluding %d voxels as active\n',nexclude);
  fprintf('Continuing with %d voxels\n',ninclude);
  y = y(:,indinclude);
  
  fprintf('Removing mean from functionals (%g)\n',toc);
  ymn = mean(y);
  y = y - repmat(ymn,[nf 1]);
  
  fprintf('Computing My (%g)\n',toc);
  My = y*y'; %'
  
  [u s blah] = svd(My);
  ds = diag(s);
  plot(100*ds/sum(ds));
  
  [dim0 pvs pvsw] = fast_estdimsvd(s);
  dim = ceil(dim0*dimfrac);
  fprintf('run = %s  dim0 = %d, dim = %d\n',runid,dim0,dim);
  
  % Only keep the first dim components of u %
  u = u(:,1:dim);

  % Test u for signal %
  beta = B*u
  uhat = X*beta;
  r = u - uhat;
  rvar = sum(r.^2)/dof;
  [F, dof1, dof2] = fast_fratio(beta,X,rvar,C);
  p = FTest(dof1, dof2, F);
  
keyboard  

end

toc









