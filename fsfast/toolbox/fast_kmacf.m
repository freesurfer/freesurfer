
% Things to do: 
%  pretaper
%  Mkjw

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

XX = load(xmatfile);
X = XX.Xfinal;
T = X*inv(X'*X)*X';
R = eye(size(T)) - T;
M = fast_kjw_mtx(R,nmaxlag);
Mkjw = inv(M);


% Fit racf with 
if(fitpolyz)
  fprintf('Fitting to polyz %d (%g)\n',polyzorder,toc);
  [pc acfmask] = fast_polyzacf_fit(acfmask,polyzorder,nmaxlag);
  fprintf(' ... done %g\n',toc);
end


yacfmasked = Mkjw*acfmask(1:nmaxlag,:);
yacfmasked = yacfmasked./repmat(yacfmasked(1,:),[nmaxlag 1]);

acfmasktrunc = yacfmasked(2:nmaxlag,:);
% apply taper to trunc only

fprintf('Starting kmeans  %g\n',toc);
[kmeans, kmap, d2min, niters, acffit] = ...
    fast_kmeans(acfmasktrunc,nclusters,[],nitersmax);
fprintf('Done kmeans  %g\n',toc);

kmapvol = zeros(1,nv);
if(~isempty(indmask))
  kmapvol(indmask) = kmap;
else
  kmapvol = kmap;
end  
kmapvol = reshape(kmapvol,[nslices nrows ncols]);
kmapstem = sprintf('%s-kmap',outbasestem);
fmri_svbvolume(kmapvol,kmapstem);

racfkm = zeros(nf,nclusters);
yacfkm = zeros(nf,nclusters);
for k = 1:nclusters
  indk = find(kmap==k);
  nindk = length(indk);
  racfk = mean(acfmask(:,indk),2); % full acf, inc 1
  yacfk = mean(yacfmasked(:,indk),2); % full acf, inc 1
  [rcnd rmineig] = fast_acfcond(racfk);
  [ycnd ymineig] = fast_acfcond(yacfk);
  fprintf('k = %2d, nindk = %5d (%4.1f%%) , ymineig = %g, ycnd=%g \n',...
	  k,nindk,100*nindk/nmask,ymineig,ycnd);
  racfkm(:,k) = racfk;
  yacfkm(1:nmaxlag,k) = yacfk;
end

tmp(1,:,:) = yacfkm';
kmeansfile = sprintf('%s-kmyacf_000.bfloat',outbasestem);
fmri_svbfile(tmp,kmeansfile);

tmp(1,:,:) = yacfkm';
kmeansfile = sprintf('%s-kmracf_000.bfloat',outbasestem);
fmri_svbfile(tmp,kmeansfile);

fprintf('fast_kmacf done %g\n',toc);

