function [prayleigh, indnz] = fast_prayleigh(vol)
% [prayleigh, indnz] = fast_prayleigh(vol)
%
% Rayleigh Probability distribution/density function sampled at
% x. If mean is not specified, it defaults to 1. The variance
% is dependent on the mean and is (4/pi-1)*mean.^2, which is
% returned in pdfvar. 
%
% Rayleigh dist noise can be constructed from complex white noise:
%   y = abs(randn(10000,1) +i*randn(10000,1));
%   The var will be the var of the white noise times (2-pi/2)
%
%


%
% fast_prayleigh.m
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

prayleigh = [];

if(nargin ~= 1)
  fprintf('[prayleigh, indnz] = fast_prayleigh(vol)\n');
  return;
end

if(isfield(vol,'vol')) 
  prayleigh = vol;
  prayleigh.vol = ones(vol.volsize);
  vol = vol.vol; 
end

vol = fast_vol2mat(vol);

volmn = mean(vol);
volstd = std(vol);
indnz = find(volstd ~= 0);
nnz = length(indnz);
indz  = find(volstd == 0);

usermse = 0;

tic
fprintf('nnz = %d\n',nnz);
for nthvox = 1:nnz
  %if(rem(nthvox,10000)==0) fprintf('nthvox = %d, %g, %g %g\n',nthvox,toc,t1,t2); end
  if(rem(nthvox,50000)==0) fprintf('nthvox = %d, %g\n',nthvox,toc); end
  v = vol(:,indnz(nthvox));
  if(~usermse)
    indok = find(v ~= 0);
    nok = length(indok);
    p = -log10(pdf_rayleigh(v(indok),mean(v)));
    pmean = mean(p);
    prayleigh.vol(indnz(nthvox)) = pmean;
  else
    [h x] = hist(v,10);
    vpdf  = h/trapz(x,h);
    vrpdf = pdf_rayleigh(x,mean(v));
    rmse  = mean((vpdf-vrpdf).^2);
    prayleigh.vol(indnz(nthvox)) = rmse;
  end
  %fprintf('%5d %g %g\n',nthvox,t1,t2);
end

prayleigh.vol(indz) = max(prayleigh.vol(indnz));

return;

  %p = pdf_rayleigh(v(indok),mean(v),1);
