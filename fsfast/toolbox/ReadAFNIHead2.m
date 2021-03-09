function headinfo =  ReadAFNIHead2(hname)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function: ReadAFNIHead2
% given the name of an AFNI HEAD
% file this function reads the 
% contents and returns information
% in a structure where the fields 
% are arranged as follows:
% contents of the structure are:
% headinfo.xdim ( x dimension size )
% headinfo.ydim ( y dimension size )
% headinfo.zdim ( z dimension size )
% headinfo.xvox ( x voxel size )
% headinfo.yvox ( y voxel size )
% headinfo.zvox ( z voxel size )
% headinfo.nvols ( number of sub-briks, i.e. for 3D+time )
% headinfo.orient ( orientation code, e.g. 351 )
% the orientation code is 3 numbers and corresponds to the following convention
% headinfo.type (type string, eg short, float)
% headinfo.byteorder (byte order string: 'b' for big endian, 'l' for little
% 0 = R-L
% 1 = L-R
% 2 = P-A
% 3 = A-P
% 4 = I-S
% 5 = S-I
% 
% Timothy M. Ellmore (tellmore@nih.gov)
% Laboratory of Brain and Cognition
% National Institute of Mental Health
%
% Modified by Douglas N. Greve
% greve@nmr.mgh.harvard.edu
% MGH-NMR Center, Boston Mass
% Added byteorder and data type


%
% ReadAFNIHead2.m
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

headinfo = [];

if(nargin ~= 1)
  fprintf('USAGE: V = ReadAFNIHead2(fname)\n');
  return;
end

fid=fopen(hname);
if(fid == -1)
  fprintf('ERROR: cannot open %s\n',hname);
  return;
end

count = 1;
while 1
 line = fgetl(fid);
 if ~isstr(line), break, end
 
 if(strcmp('name  = DELTA',line))
  voxline = count + 2;
 end

 if(strcmp('name = DATASET_DIMENSIONS',line))
   dimline = count + 2;
 end

 if(strcmp('name = BRICK_TYPES',line))
   nvolsline = count + 1;
 end

 if(strcmp('name = ORIENT_SPECIFIC',line))
   orientline = count + 2;
 end

 if(strcmp('name = BYTEORDER_STRING',line))
   byteorderline = count + 2;
 end

 % disp(line)
 count = count + 1;
end
fclose(fid);

fid=fopen(hname);
count = 1;
while 1
 line = fgetl(fid);
 if ~isstr(line), break, end
 
 if(count == voxline)
   voxinfo = sscanf(line,'%f %f %f');
 end

 if(count == dimline)
  diminfo = sscanf(line,'%d %d %d %d %d');
 end

 if(count == nvolsline)
   nvolsinfo = sscanf(line,'%s = %d');
   count = count + 1;
   line = fgetl(fid);
   datatype = sscanf(line,'%d',1);
 end

 if(count == orientline)
  orientinfo = sscanf(line,'%d %d %d');
 end
 
 if(count == byteorderline)
   byteorderstring = sscanf(line,'%s');
 end
 
 count = count + 1;
end
fclose(fid);

headinfo.xdim = diminfo(1);
headinfo.ydim = diminfo(2);
headinfo.zdim = diminfo(3);
headinfo.xvox = abs(voxinfo(1));
headinfo.yvox = abs(voxinfo(2));
headinfo.zvox = abs(voxinfo(3));
headinfo.nvols = nvolsinfo(length(nvolsinfo));

headinfo.orient = orientinfo;

headinfo.type = '';
switch datatype
  case 0, headinfo.type = 'uchar';  % Not sure about this one
  case 1, headinfo.type = 'short';
  case 3, headinfo.type = 'float';
end

headinfo.byteorder = '';
byteorderstring = byteorderstring(2:length(byteorderstring));
switch byteorderstring
  case 'LSB_FIRST~', headinfo.byteorder = 'l'; 
  case 'MSB_FIRST~', headinfo.byteorder = 'b'; 
end

return;
