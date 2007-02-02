% fast_fmriqa.m
% d = voxel size
%   gstd = d/sqrt(-4*log(ar1));
%   fwhm = gstd*sqrt(log(256.0));
%   fwhm = sqrt(log(256.0))*d/sqrt(-4*log(ar1));

tic;

m = MRIread('maske4');
indm = find(m.vol);

%infile = 'fqa.fa10.nii';
%outdir = 'qa.fa10';
infile = 'f2.nii';
outdir = 'qa';

polyorder = 2;
synth = 0;

mkdirp(outdir);

fprintf('Loading data   %g\n',toc);
f = MRIread(infile);
%f.vol = f.vol(:,:,:,1:50);
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
res = f;
res.vol = fast_mat2vol(rmat,f.volsize);

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
	mean(swnvar.vol(indm)),...
	mean(scnvar.vol(indm)));

figure(1);
imagesc(vol2mos(swnvar.vol),[-5 40]);colorbar
figure(2);
imagesc(vol2mos(scnvar.vol),[-5 5]);colorbar
figure(3);
imagesc(vol2mos(pctscnvar.vol),[0 40]);colorbar

fprintf('Done   %g\n',toc);
