% fast_fmriqa.m
% d = voxel size
%   gstd = d/sqrt(-4*log(ar1));
%   fwhm = gstd*sqrt(log(256.0));
%   fwhm = sqrt(log(256.0))*d/sqrt(-4*log(ar1));

nskip = 5;

t1  = 930; 
t2s = 47;
tr  = 2000;
te  = 30;
s90 = ssbloch(tr,te,90*pi/180,t1,t2s);   %??? 
s77 = ssbloch(tr,te,77*pi/180,t1,t2s)/s90;    
s10 = ssbloch(tr,te,10*pi/180,t1,t2s)/s90;


for flip_angle = [10 77]

infile = sprintf('f%2d.nii',flip_angle);
outdir = sprintf('qa%2d',flip_angle);
fprintf('infile %s\n',infile);

tic;

m0 = MRIread('mask');
m9 = MRIread('maske4');

indm0 = find(m0.vol);
indm9 = find(m9.vol);
nmask9 = length(indm9);


polyorder = 2;
synth = 0;

mkdirp(outdir);
fname = sprintf('%s/mask.nii',outdir);
MRIwrite(m0,fname);
fname = sprintf('%s/maske9.nii',outdir);
MRIwrite(m9,fname);


fprintf('Loading data   %g\n',toc);
f = MRIread(infile);
f.vol = f.vol(:,:,:,nskip:end);
f.nframes = size(f.vol,4);
nf = f.nframes;
nv = prod(size(f.vol));

if(synth)
  fprintf('Sythesizing %g\n',toc);
  f.vol = randn(size(f.vol)); % synth for testing
  nc = randn(1,1,1,nf);
  nc = nc/std(nc(:));
  f.vol = f.vol + sqrt(.02)*repmat(nc,[f.volsize 1]); % 2% of
else
  fprintf('NOT Sythesizing %g\n',toc);
end

X = fast_polytrendmtx(1,nf,1,polyorder);

fmat = fast_vol2mat(f);
[betamat rvarmat vdof rmat] = fast_glmfit(fmat,X);

meanmn9 = mean(betamat(1,indm9));
rvarmn9 = mean(rvarmat(indm9));
snrmn9 = meanmn9/sqrt(rvarmn9);
fprintf('mn = %g  rvar = %g  snr = %g (nvox = %d)\n',...
	meanmn9,rvarmn9,snrmn9,nmask9);

res = f;
res.vol = fast_mat2vol(rmat,f.volsize);

[U S V] = fast_svd(rmat(:,indm0));
ds2 = diag(S).^2;
pvs = 100*ds2/sum(ds2);
cpvs = cumsum(pvs);
Vtmp = zeros(size(rmat));
Vtmp(:,indm0) = V';
vv = fast_mat2vol(Vtmp,f.volsize);
fname = sprintf('%s/V.nii',outdir);
tmp = f;
tmp.vol = vv(:,:,:,1:20);
MRIwrite(tmp,fname);

fname = sprintf('%s/U.dat',outdir);
fp = fopen(fname,'w');
fmt = repmat('%g ',[1 20]);
fmt = [fmt '\n'];
fprintf(fp,fmt,U(:,1:20)');
fclose(fp);

fname = sprintf('%s/pvs.dat',outdir);
fp = fopen(fname,'w');
fprintf(fp,'%5.2f   %5.2f\n',[pvs cpvs]');
fclose(fp);

fmn = f;
fmn.vol = fast_mat2vol(betamat(1,:),f.volsize);
fname = sprintf('%s/mean.nii',outdir);
MRIwrite(fmn,fname);


fvar = f;
fvar.vol = fast_mat2vol(rvarmat,f.volsize);
fname = sprintf('%s/var.nii',outdir);
MRIwrite(fvar,fname);

fstd = fvar;
fstd.vol = sqrt(fvar.vol);
fname = sprintf('%s/std.nii',outdir);
MRIwrite(fstd,fname);

indz = find(fstd.vol < 1e-7 | fmn.vol == 0);
fsnr = fvar;
fsnr.vol = fmn.vol ./ fstd.vol;
fsnr.vol(indz) = 0;
fname = sprintf('%s/snr.nii',outdir);
MRIwrite(fsnr,fname);

fsnr2 = fsnr;
fsnr2.vol = fsnr.vol.^2;
fname = sprintf('%s/snr2.nii',outdir);
MRIwrite(fsnr2,fname);

fprintf('Normalizing residuals');
rnormmat = rmat ./ repmat(sqrt(rvarmat),[nf 1]);
rnormmat(:,indz) = 0;
rnorm = f;
rnorm.vol = fast_mat2vol(rnormmat,f.volsize);

fprintf('AR1 Row   %g\n',toc);
ar1r = f;
ar1r.vol = zeros(f.volsize);
ar1r.vol(1:end-1,:,:) = mean(rnorm.vol(1:end-1,:,:,:) .* rnorm.vol(2:end,:,:,:),4);
fname = sprintf('%s/ar1r.nii',outdir);
MRIwrite(ar1r,fname);

fprintf('AR1 Col   %g\n',toc);
ar1c = f;
ar1c.vol = zeros(f.volsize);
ar1c.vol(:,1:end-1,:) = mean(rnorm.vol(:,1:end-1,:,:) .* rnorm.vol(:,2:end,:,:),4);
fname = sprintf('%s/ar1c.nii',outdir);
MRIwrite(ar1c,fname);

fprintf('AR1 Slice   %g\n',toc);
ar1s = f;
ar1s.vol = zeros(f.volsize);
ar1s.vol(:,:,1:end-1) = mean(rnorm.vol(:,:,1:end-1,:) .* rnorm.vol(:,:,2:end,:),4);
fname = sprintf('%s/ar1s.nii',outdir);
MRIwrite(ar1s,fname);

end

return

figure(1);
imagesc(vol2mos(fvar.vol));colorbar

return

% Begin smoothing:
fprintf('Smoothing   %g\n',toc);
fwhmlist = [1 2 3]';
nfwhm = length(fwhmlist);
delta = zeros(f.volsize);
delta(round(end/2),round(end/2),round(end/2)) = 1;
fsmvar = f;
fsmvar.vol = zeros([f.volsize nfwhm]);
clear sumk2;
nth = 0;
nfuse = nf;
if(nf > 50) nfuse = 50; end
for fwhm = fwhmlist'
  fprintf('%d/%d  %g   %g\n',nth,nfwhm,fwhm,toc);
  nth = nth + 1;
  rfwhm = fwhm/f.volres(2);
  cfwhm = fwhm/f.volres(1);
  sfwhm = fwhm/f.volres(3);
  volsm = fast_smooth3d(res.vol(:,:,:,1:nfuse),cfwhm,rfwhm,sfwhm);
  fsmvar.vol(:,:,:,nth) = sum(volsm.^2,4)/(nfuse-polyorder);
  deltasm = fast_smooth3d(delta,cfwhm,rfwhm,sfwhm);
  sumk2(nth,1) = sum(deltasm(:).^2);
end

tmp = fast_vol2mat(fsmvar);
XX = [ones(nfwhm,1) sumk2];
b = inv(XX'*XX)*XX'*tmp;

scnvar = f;
scnvar.vol = fast_mat2vol(b(1,:),f.volsize);
fname = sprintf('%s/scnvar.nii',outdir);
MRIwrite(scnvar,fname);

swnvar = f;
swnvar.vol = fast_mat2vol(b(2,:),f.volsize);
fname = sprintf('%s/swnvar.nii',outdir);
MRIwrite(swnvar,fname);

pctscnvar = f;
pctscnvar.vol = 100*scnvar.vol./swnvar.vol;
pctscnvar.vol(indz) = 0;
fname = sprintf('%s/pctscnvar.nii',outdir);
MRIwrite(pctscnvar,fname);

fprintf('swn = %g, scn = %g\n', ...
	mean(swnvar.vol(indm9)),...
	mean(scnvar.vol(indm9)));

figure(1);
imagesc(vol2mos(swnvar.vol),[-5 40]);colorbar
title('SWN');
figure(2);
imagesc(vol2mos(scnvar.vol),[-5 5]);colorbar
title('SCN');
figure(3);
imagesc(vol2mos(pctscnvar.vol),[0 40]);colorbar
title('Percent SCN');

fprintf('Done   %g\n',toc);
