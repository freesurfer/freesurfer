% tdr_recon_sim.m


%
% tdr_recon_sim.m
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

tic;

rcolmat = '/home/greve/sg1/dng072203/R20.1.tik05.mat';

fprintf('Loading ... ');
load(rcolmat);
fprintf(' Done (%g)\n',toc);

nstd = 0.00005;
[nrows ncols nslices] = size(fidvol1);
nv = prod([nrows ncols nslices]);
nkcols = size(kepiref_dist,2);

Frow = fast_dftmtx(nkcols);
Frow = Frow(:,colkeep);
Rrow = inv(Frow'*Frow)*Frow';
Rrow = transpose(Rrow);

epi0 = zeros(nrows,ncols,nslices);
for slice = 1:nslices
  %fprintf('Slice %d  (%g)\n',slice,toc);
  kslice0 = kepiref_dist(:,:,slice);
  kepi02 = kslice0 * Rrow;
  for col = 1:ncols
    Rcol = Rtdr(:,:,col,slice);
    epi0(:,col,slice) = abs(Rcol * kepi02(:,col));
  end
end

nsim = 20;
epi = zeros(nrows,ncols,nslices,nsim);
for n = 1:nsim
  fprintf('n=%d, (%g)\n',n,toc);
  for slice = 1:nslices
    kslice0 = kepiref_dist(:,:,slice);
    knoise = nstd * ( randn(nrows,nkcols) + i*randn(nrows,nkcols));
    kslice = kslice0 + knoise;
    
    kepi2  = kslice * Rrow;
    for col = 1:ncols
      Rcol = Rtdr(:,:,col,slice);
      epi(:,col,slice,n)  = abs(Rcol * kepi2(:,col));
    end
  end

end

epistd = std(epi,[],4);
epimn = mean(epi,4);
mn = mean(reshape1d(epimn));
r = epistd./epimn;
ind = find(epimn < mn);
r(ind) = 0;
mos = vol2mos(r);            

figure(2)
imagesc(mos);colorbar;

return


swapview('-init','-v1',epi0,'-v2',epi(:,:,:,1));











