
% Things to do: 
%  pretaper
%  Mkjw


%
% fast_kmacf.m
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

if(0)
fsdpath = '/home/greve/sg1/dngrest/bold';
acfstem = 'faces-mcext-004/acf/acf001';
maskstem = 'masks/brain';
%outstem = 'faces-mcext-004/acf/acf001';
outbasestem = '';
nitersmax = 10;
nmaxlag = 20;
nclusters = 10;
cd(fsdpath);
end

tic;

if(isempty(outbasestem))
  outbasestem = acfstem;
end


fprintf('nclusters = %d\n',nclusters);
fprintf('nmaxlag   = %d\n',nmaxlag);
fprintf('nitersmax = %d\n',nitersmax);

fprintf('Loading %s  %g\n',acfstem,toc);
acf = fmri_ldbvolume(acfstem);
fprintf('Done loading %s  %g\n',acfstem,toc);
if(isempty(acf)) return; end
[nslices nrows ncols nf] = size(acf);
nv = nslices*nrows*ncols;
acf = reshape(acf,[nv nf])';

if(~isempty(maskstem))
  mask = fmri_ldbvolume(maskstem);
  if(isempty(mask)) return; end
  indmask = find(mask);
  nmask = length(indmask);
  if(isempty(indmask))
    fprintf('ERROR: mask %s has no voxels\n');
    return;
  end
  fprintf('nmask = %d\n',nmask);
  acfmask = acf(:,indmask);
else
  indmask = [];
  nmask = nv;
  acfmask = acf;
end

if(nmaxlag == -1 | nmaxlag > nf) nmaxlag = nf; end

% Fit racf with polyz
if(fitpolyz)
  fprintf('Fitting to polyz %d (%g)\n',polyzorder,toc);
  [pc acfmask] = fast_polyzacf_fit(acfmask,polyzorder,nmaxlag);
  fprintf(' ... done %g\n',toc);
end

if(usekjw)
  XX = load(xmatfile);
  X = XX.Xfinal;
  T = X*inv(X'*X)*X';
  R = eye(size(T)) - T;
  M = fast_kjw_mtx(R,nmaxlag);
  Mkjw = inv(M);
  yacfmasked = Mkjw*acfmask(1:nmaxlag,:);
  yacfmasked = yacfmasked./repmat(yacfmasked(1,:),[nmaxlag 1]);
  acfmasktrunc = yacfmasked(2:nmaxlag,:);
else
  acfmasktrunc = acfmask(2:nmaxlag,:);
end

kmeans0 = []; 
if(~isempty(initmethod) & strcmpi(initmethod,'firstlag'))
  fprintf('Initializing kmeans with first lag %g\n',toc);
  [h x] = hist(acfmasktrunc(1,:),nclusters);
  [kmeans, kmap, d2min, niters, acffit] = ...
      fast_kmeans(acfmasktrunc(1,:),nclusters,x,10);
  kmeans0 = zeros(nmaxlag-1,nclusters);
  for k = 1:nclusters
    indk = find(kmap==k);
    nindk = length(indk);
    %fprintf('k=%d   %d\n',k,nindk);
    racfk = mean(acfmasktrunc(:,indk),2);
    kmeans0(:,k) = racfk;
  end
end
if(~isempty(initmethod) & strcmpi(initmethod,'svd'))
  fprintf('Initializing kmeans with svd %g\n',toc);
  Macf = acfmask *acfmask';
  [u s v] = svd(Macf);
  kmeans0 = u(1:nmaxlag-1,1:nclusters);
end

fprintf('Starting kmeans  %g\n',toc);
[kmeans, kmap, d2min, niters, acffit] = ...
    fast_kmeans(acfmasktrunc,nclusters,kmeans0,nitersmax);
fprintf('Done kmeans  %g\n',toc);

% Sort so that the racf closest to white is first %
kmeans0 = kmeans;
kmap0 = kmap;
rms = sum(kmeans.^2);
[rmssorted indsort] = sort(rms);
kmeans = kmeans(:,indsort);
kmaptmp = kmap0;
for k = 1:nclusters
  indk = find(kmap==indsort(k));  
  kmaptmp(indk) = k;
end
kmap = kmaptmp;

racfkm = zeros(nf,nclusters);
if(usekjw) yacfkm = zeros(nf,nclusters); end
for k = 1:nclusters
  indk = find(kmap==k);
  nindk = length(indk);
  racfk = mean(acfmask(:,indk),2); % full acf, inc 1
  [rcnd rmineig] = fast_acfcond(racfk);
  racfkm(:,k) = racfk;
  fprintf('k = %2d, nindk = %5d (%4.1f%%) %6.5f ',...
	  k,nindk,100*nindk/nmask,sqrt(mean(racfk(2:end).^2)));
  if(usekjw) 
    yacfk = mean(yacfmasked(:,indk),2);  % full acf, inc 1
    [ycnd ymineig] = fast_acfcond(yacfk); 
    yacfkm(1:nmaxlag,k) = yacfk;
    fprintf(' ymineig = %g, ycnd=%g \n',ymineig,ycnd);
  else
    fprintf('\n');
  end
end

kmapvol = zeros(1,nv);
if(~isempty(indmask)) kmapvol(indmask) = kmap;
else                  kmapvol = kmap;
end  
kmapvol = reshape(kmapvol,[nslices nrows ncols]);
kmapstem = sprintf('%s-%s-kmap',outbasestem,kmacfid);
fprintf('Saving kmap to %s\n',kmapstem);
fmri_svbvolume(kmapvol,kmapstem);

clear tmp;
tmp(1,:,:) = racfkm';
kmeansfile = sprintf('%s-%s-racf_000.bfloat',outbasestem,kmacfid);
fprintf('Saving kmeans to %s\n',kmeansfile);
fmri_svbfile(tmp,kmeansfile);

if(usekjw)
  clear tmp;
  tmp(1,:,:) = yacfkm';
  kmeansfile = sprintf('%s-kmyacf_000.bfloat',outbasestem);
  fmri_svbfile(tmp,kmeansfile);
end

fprintf('fast_kmacf done %g\n',toc);

