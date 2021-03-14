function err = fmri_svbfile(y, BFileName, Endian)
%
% fmri_svbfile(y,BFileName,Endian)
%
% Saves a bshort or bfloat given the full path
% name of the BFile.  The type (bshort or
% bfloat) is determined from the name.
% The header is written with the dimensions.
% Converts from matlab column-major format
% to row major.
%
% Endian is either 0 or 1 (default is 0)
%
% See also: LdBFile
%
%


%
% fmri_svbfile.m
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

err = 1 ;
if(nargin ~= 2 & nargin ~= 3) 
  error('USAGE: SvBFile(y,BFileName,<Endian>)');
end

if(nargin == 2) Endian = 0; end

if(Endian ~= 0 & Endian ~= 1)
  msg = sprintf('Endian = %d, must me either 0 or 1',Endian);
  qoe(msg); error(msg);
end

BFileName = deblank(BFileName);
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
  msg = 'BFileName must have a non-null base';
  qoe(msg);
  error(msg);
end

HdrFile   = strcat(Base,'.hdr');

%%% Open the header file %%%%
fid=fopen(HdrFile,'w');
if fid == -1 
  msg = sprintf('Could not open header %s for writing\n',HdrFile);
  qoe(msg);
  error(msg);
end

ndy = length(size(y));
nR = size(y,1);
nC = size(y,2);
nD = prod(size(y))/(nR*nC);

%%%% Write the Dimension to the header %%%%
fprintf(fid,'%d %d %d %d\n',nR,nC,nD,Endian); % 0=big-endian
fclose(fid);

%%% Open the bfile %%%%
if(Endian == 0) EndianFlag = 'b';
else            EndianFlag = 'l';
end
fid=fopen(BFileName,'w',EndianFlag);
if fid == -1 
  msg = sprintf('Could not open bfile %s for writing\n',BFileName);
  qoe(msg);
  error(msg);
end

%%%% Transpose into row-major format %%%%
y = reshape(y, [nR nC nD]);
y = permute(y, [2 1 3]);

%%%%% Save the Slice %%%%%
count = fwrite(fid,y,precision);
fclose(fid); 

if(count ~= prod(size(y)))
  msg = sprintf('Did not write correct number of items (%d/%d)',...
                count,prod(size(y)));
  qoe(msg);  error(msg);
end

err = 0;

return;
