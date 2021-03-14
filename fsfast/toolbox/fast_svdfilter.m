% fast_svdfilter.m

% extregname
% mergeextreg
% mergedextreg
% filter option


%
% fast_svdfilter.m
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

tic

fsd = '/space/greve/1/users/greve/wagner/JL/bold';
analysis = 'onoff';
runidlist = strvcat('006','007','008','009');

fsd = '/home/greve/sg1/nh-rest/bold';
analysis = 'onoff60s';
runidlist = strvcat('006','007','008','009');

funcrstem = 'fmcsm5';
maskrstem = 'brain';
dimfrac   = 1;
sigthresh = 3;


%--------------------------------------------------%
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
  nf = size(X,1);
  Ntask = XX.Nnnc*XX.Navgs_per_cond;
  Nbeta = size(X,2);
  C = zeros(Ntask,Nbeta);
  C(1:Ntask,1:Ntask) = 1;
  dof = size(X,1)-size(X,2);
  B = inv(X'*X)*X';
  R = eye(nf)-X*B;
  
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
  y0 = y;
  
  if(0)
    fprintf('Computing whitening matrix (%g)\n',toc);
    r = R*y;
    racf = fast_acorr(r);
    racfmn0 = mean(racf,2);
    [cnd mineig S] = fast_acfcond(racfmn0);
    racfmn = racfmn0;
    while(mineig < 0 & cnd > 100)
      racfmn = racfmn .* tukeytaper(nf);
      [cnd mineig S] = fast_acfcond(racfmn);
    end
    [u s v] = svd(S);
    W = u*inv(sqrt(s))*v';
    WX = W*X;
    Wy = W*y;
    WB = inv(WX'*WX)*WX';
    beta = WB*Wy;
    r = Wy-WX*beta;
  else
    beta = B*y;
    r = y-X*beta;
  end
  
  fprintf('Analyzing functionals to get exclusion mask (%g)\n',toc);
  rvar = sum(r.^2)/dof;
  [F, dof1, dof2, g] = fast_fratio(beta,X,rvar,C);
  p = FTest(dof1, dof2, F);
  [pdf alpha nx fpr] = ComputePDF(p,0,1,.001);
  
  indexclude = find(p < pthresh);
  indinclude = find(p >= pthresh);
  nexclude = length(indexclude);
  ninclude = length(indinclude);
  fprintf('Excluding %d (%g%%) voxels as active (%g)\n',...
	  nexclude,100*nexclude/nmask,sigthresh);
  fprintf('Continuing with %d voxels\n',ninclude);
  y = y(:,indinclude);
  
  fprintf('Removing mean from functionals (%g)\n',toc);
  ymn = mean(y);
  y = y - repmat(ymn,[nf 1]);
  
  fprintf('Computing My (%g)\n',toc);
  My = y*y'; %'
  
  [u s blah] = svd(My);
  ds = diag(s);
  %plot(100*ds/sum(ds));
  
  [dim0 pvs pvsw] = fast_estdimsvd(s);
  dim = ceil(dim0*dimfrac);
  fprintf('run = %s  dim0 = %d (%g), dim = %d (th=%g)\n',...
	  runid,dim0,sum(pvs(1:dim0)),dim,sum(pvs(1:dim)));
  
  % Only keep the first dim components of u %
  u = u(:,1:dim);

  % Test u for signal %
  beta = B*u;
  uhat = X*beta;
  r = u - uhat;
  rvar = sum(r.^2)/dof;
  [F, dof1, dof2] = fast_fratio(beta,X,rvar,C);
  p = FTest(dof1, dof2, F);
  
  X = [X u];
  Ntask = XX.Nnnc*XX.Navgs_per_cond;
  Nbeta = size(X,2);
  C = zeros(Ntask,Nbeta);
  C(1:Ntask,1:Ntask) = 1;
  dof = size(X,1)-size(X,2);
  B = inv(X'*X)*X';
  beta = B*y0;
  r = y0-X*beta;
  rvar = sum(r.^2)/dof;
  [F, dof1, dof2, g] = fast_fratio(beta,X,rvar,C);
  p = FTest(dof1, dof2, F);
  [pdf2 alpha2 nx fpr2] = ComputePDF(p,0,1,.001);
  
  figure(nthrun);
  loglog(alpha,alpha,alpha,fpr,alpha2,fpr2);
  legend('ideal','none','svdfil');
  title(sprintf('Run %s',runid));
  
end

toc









