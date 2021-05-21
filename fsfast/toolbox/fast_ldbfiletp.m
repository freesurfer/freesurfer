function y = fast_ldbfile(bfile,frame)
%
% y = fast_ldbfile(bfile,frame)
%
% Loads a single frame from a bshort or bfloat 
% given the full path and name of the BFile.  The type 
% (bshort or bfloat is determined from the name).
%
%
%
% See also: fmri_ldbfile()


%
% fast_ldbfiletp.m
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

if(nargin ~= 1 & nargin ~= 2)
  fprintf(2,'USAGE: y = fast_ldbfile(bfile,frame)');
  qoe; return;
end

if(nargin == 1) getall = 1;
else            getall = 0;
end

ks = findstr(bfile,'.bshort');
kf = findstr(bfile,'.bfloat');

if(isempty(ks) & isempty(kf))
  msg = 'bfile must be either bshort or bfloat';
  qoe(msg); error(msg);
end

if( ~isempty(ks) ) 
  precision = 'int16';
  Base = bfile(1:ks-1);
else               
  precision = 'float32';
  Base = bfile(1:kf-1);
end

if( isempty(Base) )
  s = 'fast_ldbfile: bfile must have a non-null base';
  qoe(msg);
  error(msg);
end

%%% Open the header file %%%%
HdrFile = strcat(Base,'.hdr');
fid=fopen(HdrFile,'r');
if fid == -1 
  msg = sprintf('fast_ldbfile: Could not open %s file',HdrFile); 
  qoe(msg); error(msg);
end

%%%% Read the Dimension from the header %%%%
hdr=fscanf(fid,'%d',[1,4]);
fclose(fid);
Nr  = hdr(1);
Nc  = hdr(2);
Nf  = hdr(3);
Endian = hdr(4);

if(~getall)
  if(frame >= Nf)
    msg = printf('ERROR: frame (%d) >= nframes (%d)\n',frame,Nf);
    qoe(msg); error(msg);
  end
end

%%%% Open the bfile %%%%%
if(Endian == 0) fid=fopen(bfile,'r','b'); % Big-Endian
else            fid=fopen(bfile,'r','l'); % Little-Endian
end
if( fid == -1)
  msg = sprintf('ERROR: fast_ldbfile: Could not open %s file',bfile); 
  qoe(msg); error(msg);
end

if(getall)

  %%% Read the entire bfile %%%
  [y count] = fread(fid,[Nr*Nc*Nf],precision);
  fclose(fid); 
  if(count ~= Nr*Nc*Nf)
    msg = sprintf('Read %d from %s, expecting %d\n',...
                  count,bfile,Nr*Nc*Nf);
    qoe(msg); error(msg);  
  end

  %% Reshape into image dimensions %%
  y = reshape(y,[Nc Nr Nf]);

  %%% Transpose because matlab uses row-major %%%
  y = permute(y,[2 1 3]);

else

  %%% Seek to the desired frame %%%
  status = fseek(fid,Nr*Nc*frame,'bof');
  if(status)
    msg = sprintf('ERROR: fast_ldbfile: seeking to %dth frame\n',frame);
    qoe(msg); error(msg);  
  end

  %%% Read the single frame %%%
  [y count] = fread(fid,[Nr*Nc],precision);
  fclose(fid); 
  if(count ~= Nr*Nc)
    msg = sprintf('ERROR: fast_ldbfile: Read %d from %s, expecting %d\n',...
                  count,bfile,Nr*Nc);
    qoe(msg); error(msg);  
  end

  %% Reshape into image dimensions %%
  y = reshape(y,[Nc Nr]);

  %%% Transpose because matlab uses row-major %%%
  y = permute(y,[2 1]);

end

return;



