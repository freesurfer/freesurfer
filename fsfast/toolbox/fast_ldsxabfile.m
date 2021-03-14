function [havg, eresvar, sxadat] = fast_ldsxabfile(sxabfile)
%
% [havg eresvar sxadat] = fast_ldsxabfile(sxabfile)
%
% This function reads in the selxavg values from the given bfile
% assuming that the data are stored in selavg format.
%
% havg - (nrows,ncols,Nhtot)
% eresvar - (nrows,ncols) - residual error variance
% sxadat - info from the .dat file
%
%
%
% See also: fmri_svsxabvol()


%
% fast_ldsxabfile.m
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
  msg = 'USAGE: [havg eresvar sxadat] = fast_ldsxabfile(sxabfile)';
  qoe(msg); error(msg);
end

sxastem = fmri_getstem(sxabfile);

datfile = sprintf('%s.dat',sxastem);
sxadat = fmri_lddat3(datfile);
if(isempty(sxadat)) return; end

[nrows ncols ntp fs ns endian bext] = fmri_bfiledim(sxastem);

eresvar = zeros(nrows,ncols);
havg = zeros(nrows,ncols,sxadat.Nh*sxadat.Nnnc);
[avglist stdlist] = fast_hlist(sxadat.Nc,sxadat.Nh);

%%%% Open the bfile %%%%%
if(endian == 0) fid=fopen(sxabfile,'r','b'); % Big-Endian
else            fid=fopen(sxabfile,'r','l'); % Little-Endian
end
if fid == -1 
  msg = sprintf('LdBFile: Could not open %s file',sxabfile); 
  qoe(msg);
  error(msg);
end
precision = 'float32';

% Get past the first Nh planes to the variance plane %
nskip = 4*nrows*ncols*sxadat.Nh;
status = fseek(fid,nskip,'bof');

% Read in the variance plane %
[z count] = fread(fid,[ncols nrows],precision);
if(count ~= nrows*ncols)
  msg = sprintf('Read %d from %s, expecting %d\n',count,sxabfile,nrows*ncols);
  fclose(fid); qoe(msg); error(msg);  
end
eresvar = z'; %'
eresvar = eresvar.*eresvar;

% Get past the next (Nh-1) planes to the first average plane %
nskip = 4*nrows*ncols*(sxadat.Nh-1);
status = fseek(fid,nskip,'cof');

% Go through each non-null condition and read in the averages %
nskip = 4*nrows*ncols;
n = 1;
for c = 1:sxadat.Nnnc
  for s = 1:2
    for h = 1:sxadat.Nh
      if(s==2) status = fseek(fid,nskip,'cof');
      else
        [z count] = fread(fid,[ncols nrows],precision);
        if(count ~= nrows*ncols)
          msg = sprintf('Read %d from %s, expecting %d\n',...
                         count,sxabfile,nrows*ncols);
          fclose(fid); qoe(msg); error(msg);  
        end
        havg(:,:,n) = z'; %'
        n = n + 1;
      end
      %fprintf('%2d %d %2d\n',c,s,h);
    end
  end
end

fclose(fid);

return;


