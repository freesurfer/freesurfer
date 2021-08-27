function [havg, eresvar, sxadat] = fast_ldsxavol(sxafile)
%
% [havg eresvar sxadat] = fast_ldsxavol(sxafile)
%
% This function reads in the selxavg values from the given file
% assuming that the data are stored in selxavg format. The sxafile
% can be a stem or have .bhdr, .mgh, .mgz, .nii, or .ngz extension.
% Uses MRIread().
%
% havg - mri struct with (nrows,ncols,nslice,Nhtot)
% eresvar -  residual error variance mri struct with (nrows,ncols,nslices,1)
% sxadat - info from the .dat file
%
% The selxavg format is:
% First Nh frames are 0
% Next  Nh frames are resstd
% Next  Nh frames are regressors for first condition
% Next  Nh frames are std of reg for first condition
% Next  Nh frames are regressors for second condition
% Next  Nh frames are std of reg for second condition
% ...
%
% See also: fast_ldsxabfile(), fmri_svsxabvol(), fmri_hlist().
%
%


%
% fast_ldsxavol.m
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

havg = [];
eresvar = [];
sxadat = [];

if(nargin ~= 1) 
  msg = '[havg eresvar sxadat] = fast_ldsxabfile(sxafile)';
  qoe(msg); error(msg);
end

[fspec sxastem fmt] = MRIfspec(sxafile);
if(isempty(fspec))
  fprintf('ERROR: cannot determine format of %s\n',sxafile);
  return;
end

datfile = sprintf('%s.dat',sxastem);
sxadat = fmri_lddat3(datfile);
if(isempty(sxadat)) return; end

sxa = MRIread(fspec);
if(isempty(sxa)) return; end

eresvar = sxa;
eresvar.vol = (sxa.vol(:,:,:,sxadat.Nh+1)).^2;

% List of indices of the regression coefficients (including 0s at first)
regind = fast_hlist(sxadat.Nc,sxadat.Nh);
regind = regind(sxadat.Nh+1:end);
havg = sxa;
havg.vol = sxa.vol(:,:,:,regind);


return;


