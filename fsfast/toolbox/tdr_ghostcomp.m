% tdr_ghostcomp.m
% $Id: tdr_ghostcomp.m,v 1.1 2003/11/25 21:05:52 greve Exp $

mghdir = '/autofs/space/greve_002/users/greve/fb-104.2/sm1/mgh';
measasc = sprintf('%s/meas.asc',mghdir);
thresh = 2;
order = 2;
rvarthresh = 0.25;
outmatfile = sprintf('%s/ghostcomp.mat',mghdir); 

% Read in the epi timing parameters
echospacing = tdr_measasc(measasc,'m_lEchoSpacing');
tDwell = tdr_measasc(measasc,'sRXSPEC.alDwellTime[0]')/1000; %usec
tRampUp = tdr_measasc(measasc,'m_alRegridRampupTime');
tFlat   = tdr_measasc(measasc,'m_alRegridFlattopTime');
tRampDown = tdr_measasc(measasc,'m_alRegridRampdownTime');
tDelSamp = tdr_measasc(measasc,'m_alRegridDelaySamplesTime');

kpcnrfile = sprintf('%s/pcnr.mgh',mghdir);
kpcnr = load_mgh(kpcnrfile);
if(isempty(kpcnr)) return; end

kpcnifile = sprintf('%s/pcni.mgh',mghdir);
kpcni = load_mgh(kpcnifile);
if(isempty(kpcni)) return; end

kpcn = kpcnr + i*kpcni;
kpcn = permute(kpcn,[2 1 3 4]);
kpcn(2,:,:,:) = flipdim(kpcn(2,:,:,:),2);

[npcn ncols nslices nframes] = size(kpcn);

[kvec gvec] = kspacevector2(ncols,tDwell,tRampUp,tFlat,...
                              tRampDown,tDelSamp,0);

% Compute the Ideal col and row DFT reconstruction matrices
Frow = fast_dftmtx(kvec);
Frow = fast_svdregpct(Frow,90);
Rrow = transpose(inv(Frow));

kpcntmp = transpose(squeeze(kpcn(1,:,:,1)));
anat = abs(kpcntmp*Rrow);

beta   = zeros(order+1,nslices,nframes);
resvar = zeros(nslices,nframes);
nover  = zeros(nslices,nframes);
thetahat = zeros(nslices,ncols,nframes);
for sliceno = 1:nslices

  kpcn1 = squeeze(kpcn(1,:,sliceno,:));
  kpcn2 = squeeze(kpcn(2,:,sliceno,:));
  kpcn3 = squeeze(kpcn(3,:,sliceno,:));
  kpcn1 = (kpcn1+kpcn3)/2;
  pcn1 = (kpcn1.'*Rrow).';
  pcn2 = (kpcn2.'*Rrow).';
    
  [thetahattmp betatmp rvartmp novertmp] = ...
      tdr_phdist_est(pcn1,pcn2,thresh,order);
  beta(:,sliceno,:) = betatmp;
  resvar(sliceno,:) = rvartmp;
  nover(sliceno,:) = novertmp;
  thetahat(sliceno,:,:) = thetahattmp;

  if(0)
    betaslice = squeeze(beta(:,sliceno,:)); 
    fprintf('slice %2d: mean  ',sliceno);
    fprintf('%9.4f ',mean(betaslice,2));
    %fprintf('\n');
    fprintf('   std  ');
    fprintf('%9.4f ',std(betaslice,[],2))
    fprintf('\n');
  end
  
end

resvarmn = mean(resvar,2);
resvargmn = mean(resvarmn);
indslices = find(resvarmn < resvargmn*rvarthresh);
thetahatmn = mean(mean(thetahat(indslices,:,:),1),3);

%-----------------------------------------------%
% This just allows a check %
kpcn1 = transpose(squeeze(kpcn(1,:,:,1)));
kpcn2 = transpose(squeeze(kpcn(2,:,:,1)));
pcn2 = kpcn2*Rrow;
pcn1hat = pcn2 .* repmat(exp(i*thetahatmn),[nslices 1]);
kpcn1hat = pcn1hat*Frow;

save(outmatfile);

%return;
%-------------------------------------------------------%
%-------------------------------------------------------%
%-------------------------------------------------------%

nn = 1:ncols;
sliceno = 1;
plot(nn,abs(kpcn1(sliceno,:)),nn,abs(kpcn1hat(sliceno,:)))

%-------------------------------------------------------%
kepirfile = sprintf('%s/echo001r.mgh',mghdir);
kepiifile = sprintf('%s/echo001i.mgh',mghdir);

vol1 = zeros(64,64,35);
vol2 = zeros(64,64,35);
vol3 = zeros(64,64,35);

for sliceno = 1:35
kepir = load_mgh(kepirfile,sliceno,1);
if(isempty(kepir)) return; end

kepii = load_mgh(kepiifile,sliceno,1);
if(isempty(kepii)) return; end

kepi = kepir + i*kepii;
kepi = permute(kepi,[2 1 3 4]);
[nrows ncols nslices nframes] = size(kepi);
oddrows  = [1:2:nrows];
noddrows = length(oddrows);
evenrows = [2:2:nrows];
nevenrows = length(evenrows);

kepi(evenrows,:,:,:) = flipdim(kepi(evenrows,:,:,:),2);
%kepi(evenrows,:) = fliplr(kepi(evenrows,:));

Rcol = inv(fast_dftmtx(nrows));

ph = thetahatmn;
ph = ph - mean(ph);
ph = ph + mean(thetahat(sliceno,:,1));
epih = kepi * Rrow;
epih(evenrows,:,:,:) = epih(evenrows,:,:,:) .* ...
    repmat(exp(i*ph/2),[nevenrows 1 nslices nframes]);
epih(oddrows,:,:,:) = epih(oddrows,:,:,:) .* ...
    repmat(exp(-i*ph/2),[noddrows 1 nslices nframes]);
kepihat = epih*Frow;
epi1 = abs(Rcol*epih);
epi1 = epi1(:,32:32+63);
vol1(:,:,sliceno) = epi1;
%figure(1);imagesc(epi1);colorbar;axis image
%colormap gray;
%cmax = max(reshape1d(epi1));
%caxis([0 cmax]);

ph = thetahat(sliceno,:,1);
epih = kepi * Rrow;
epih(evenrows,:,:,:) = epih(evenrows,:,:,:) .* ...
    repmat(exp(i*ph/2),[nevenrows 1 nslices nframes]);
epih(oddrows,:,:,:) = epih(oddrows,:,:,:) .* ...
    repmat(exp(-i*ph/2),[noddrows 1 nslices nframes]);
kepihat = epih*Frow;
epi2 = abs(Rcol*epih);
epi2 = epi2(:,32:32+63);
vol2(:,:,sliceno) = epi2;
%figure(2);imagesc(epi2);colorbar;axis image
%colormap gray;
%caxis([0 cmax]);


[kimg2, betatmptmp] = tdr_deghost(kepi,Rcol,Rrow,0);
epi3 = abs(Rcol*kimg2);
epi3 = epi3(:,32:32+63);
vol3(:,:,sliceno) = epi3;
%figure(3);imagesc(epi3);colorbar;axis image
%colormap gray;
%caxis([0 cmax]);

end
