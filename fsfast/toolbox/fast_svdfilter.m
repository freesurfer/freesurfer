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

instem = '/home/greve/sg1/dngrest/bold/004/f';
maskstem = '/home/greve/sg1/dngrest/bold/masks/brain';
xmatfile = '/home/greve/sg1/dngrest/bold/onoff/X004.mat';
pforder = -1;
pthresh = .01;
pvsthresh = 0;

XX = load(xmatfile);
if(isempty(XX))
  fprintf('ERROR: could not load %s\n',xmatfile);
  return; quit;
end
X = XX.Xfinal;

mask = fmri_ldbvolume(maskstem);
if(isempty(mask))
  fprintf('ERROR: could not load %s\n',maskstem);
  return; quit;
end
indmask = find(mask);
nmask = length(find(mask));
if(nmask == 0)
  fprintf('ERROR: could not find any voxels in mask\n');
  return; quit;
end

y = fmri_ldbvolume(instem);
if(isempty(y))
  fprintf('ERROR: could not load %s\n',instem);
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

fprintf('INFO: input vol size %d %d %d %d\n',size(y));
fprintf('INFO: found %d voxels in mask\n',nmask);

y = reshape(y,[nv nf])'; %'
y = y(:,indmask);

if(1)
  Xpf = fast_polytrendmtx(1,nf,1,XX.pfOrder);
  Epf = eye(nf) - Xpf*inv(Xpf'*Xpf)*Xpf';
  y= Epf*y;
end

fprintf('Computing My\n');
My = y*y'; %'

[u s blah] = svd(My);
ds = diag(s);
plot(100*ds/sum(ds));


dim = fast_estdimsvd(s);



toc









