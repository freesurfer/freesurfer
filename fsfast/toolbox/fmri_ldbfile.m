function y = fmri_ldbfile(varargin)
%
% y = fmri_ldbfile(bfilename)
% y = fmri_ldbfile(bfilename1,bfilename2,...,bfilenameN)
%
% Loads a bshort or bfloat given the full path
% and name of the BFile.  The type (bshort or
% bfloat is determined from the name).
% The header is read to get the dimensions, and
% the image, y, is reshaped so that it is of the correct
% dimensionality. Converts from row-major to matlabs 
% column-major.  If multiple bfiles are specified, then
% another dimension is added to y at the end to indicate
% the file from which it came.  Data from all files must
% have the same dimensionality.
%
%
%
% See also: fmri_svbile()


%
% fmri_ldbfile.m
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

y = [];

if(nargin == 0) 
  fprintf(2,'USAGE: LdBFile(BFileName)');
  qoe;
  return;
end

if( length(varargin) == 1)
  BFileList = varargin{1};
  nRuns = size(BFileList,1);
else
  nRuns = length(varargin);
  BFileList = '';
  for r = 1:nRuns,
    BFileList = strvcat(BFileList,varargin{r});
  end
end

for r = 1:nRuns,

  BFileName = deblank(BFileList(r,:));
  ks = findstr(BFileName,'.bshort');
  kf = findstr(BFileName,'.bfloat');

  if(isempty(ks) & isempty(kf))
    msg = 'BFileName must be either bshort or bfloat';
    qoe(msg);
    error(msg);
  end

  if( ~isempty(ks) ) 
    precision = 'int16';
    Base = BFileName(1:ks-1);
  else               
    precision = 'float32';
    Base = BFileName(1:kf-1);
  end

  if( isempty(Base) )
    s = 'LdBFile: BFileName must have a non-null base';
    qoe(msg);
    error(msg);
  end

  %%% Open the header file %%%%
  HdrFile = strcat(Base,'.hdr');
  fid=fopen(HdrFile,'r');
  if fid == -1 
    msg = sprintf('LdBFile: Could not open %s file',HdrFile); 
    qoe(msg);
    error(msg);
  end

  %%%% Read the Dimension from the header %%%%
  hdr=fscanf(fid,'%d',[1,4]);
  fclose(fid);
  nR  = hdr(1);
  nC  = hdr(2);
  nD  = hdr(3);
  Endian = hdr(4);

  %%%% Open the bfile %%%%%
  if(Endian == 0) fid=fopen(BFileName,'r','b'); % Big-Endian
  else            fid=fopen(BFileName,'r','l'); % Little-Endian
  end
  if fid == -1 
    msg = sprintf('LdBFile: Could not open %s file',BFileName); 
    qoe(msg);
    error(msg);
  end

  %%% Read the file in bfile %%%
  [z count] = fread(fid,[nR*nC*nD],precision);
  fclose(fid); 
  if(count ~= nR*nC*nD)
    msg = sprintf('Read %d from %s, expecting %d\n',...
                  count,BFileName,nR*nC*nD);
    qoe(msg); error(msg);  
  end

  %% Reshape into image dimensions %%
  z = reshape(z,[nC nR nD]);

  %%% Transpose because matlab uses row-major %%%
  z = permute(z,[2 1 3]);

  if(size(z,1) == 1 & size(z,3) == 1)
    y(:,:,r) = z;
  else
    y(:,:,:,r) = z;
  end

end

return;

%%% y now has size(y) = [nR nC nD nRuns] %%%

