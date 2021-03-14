function [T, VolRes, VolDim, ci] = load_corinfo(corinfofile);
%
% [T VolRes VolDim ci] = load_corinfo(corinfofile);
%
% Parses the COR-.info header file for COR volumes.
%
% T is the affine transform that converts col, row, slice to
% x, y, z in scanner coordinates, assuming that the count starts
% at 0, ie,
%
%      x         col
%      y  = T *  row
%      z        slice
%      1          1
%
% VolRes is a 3x1 vector with the col resolution, row
% resolution, and slice resolution
%
% VolDim is a 3x1 vector with the number of columns, rows, and slices.
%
% ci is a structure whose components are named after that found in 
% the COR-.info file.
%


%
% load_corinfo.m
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

if(nargin ~= 1)
  msg = '[T VolRes VolDim ci] = load_corinfo(corinfofile)';
  qoe(msg); error(msg);
end

fid = fopen(corinfofile,'r');
if(fid == -1)
  msg = sprintf('Could not open %s',corinfofile);
  qoe(msg); error(msg);
end

line = fgetl(fid); ci.imnr0 = sscanf(line,'%*s %d',1);
line = fgetl(fid); ci.imnr1 = sscanf(line,'%*s %d',1);
line = fgetl(fid); ci.ptype = sscanf(line,'%*s %d',1);
line = fgetl(fid); ci.x     = sscanf(line,'%*s %d',1);
line = fgetl(fid); ci.y     = sscanf(line,'%*s %d',1);
line = fgetl(fid); ci.fov   = sscanf(line,'%*s %d',1);
line = fgetl(fid); ci.thick = sscanf(line,'%*s %f',1);
line = fgetl(fid); ci.psiz  = sscanf(line,'%*s %f',1);
line = fgetl(fid); ci.locatn = sscanf(line,'%*s %f',1);
line = fgetl(fid); ci.strtx = sscanf(line,'%*s %f',1);
line = fgetl(fid); ci.endx  = sscanf(line,'%*s %f',1);
line = fgetl(fid); ci.strty = sscanf(line,'%*s %f',1);
line = fgetl(fid); ci.endy  = sscanf(line,'%*s %f',1);
line = fgetl(fid); ci.strtz = sscanf(line,'%*s %f',1);
line = fgetl(fid); ci.endz  = sscanf(line,'%*s %f',1);
line = fgetl(fid); ci.tr    = sscanf(line,'%*s %f',1);
line = fgetl(fid); ci.te    = sscanf(line,'%*s %f',1);
line = fgetl(fid); ci.ti    = sscanf(line,'%*s %f',1);
line = fgetl(fid); ci.flip  = sscanf(line,'%*s %f',1);

line = fgetl(fid); 
tag = sscanf(line,'%s',1);
ci.xform = [];
if(strmatch(tag,'xform'))
  ci.xform = sscanf(line,'%*s %s',1);
  line = fgetl(fid); 
end

ci.ras_good_flag = sscanf(line,'%*s %f',1);

line = fgetl(fid); ci.x_ras  = sscanf(line,'%*s %f %f %f',3);
line = fgetl(fid); ci.y_ras  = sscanf(line,'%*s %f %f %f',3);
line = fgetl(fid); ci.z_ras  = sscanf(line,'%*s %f %f %f',3);
line = fgetl(fid); ci.c_ras  = sscanf(line,'%*s %f %f %f',3);

% Direction Cosines for Col, Row, and Slice
Vc = ci.x_ras;
Vr = ci.y_ras;
Vs = ci.z_ras;

col_res = 1000*ci.psiz;
row_res = 1000*ci.psiz;
slice_res = 1000*ci.thick;
VolRes = [col_res row_res slice_res]'; %'

ncols = ci.x;
nrows = ci.y;
nslices = ci.imnr1 - ci.imnr0 + 1;
VolDim = [ncols nrows nslices]'; %'

% Construct the matrix of Direction Cosines %
Mdc = [Vc Vr Vs];

% xyz of the center of the volume
Pxyz = ci.c_ras;

% crs of the center of the volume
Pcrs = ([ncols nrows nslices]'-1)/2; %'

% Compute the xyz at crs = 0 
D3 = diag([col_res row_res slice_res]);
Pxyz0 = Pxyz - Mdc*D3*Pcrs;

% Compute the final transform
D = diag([col_res row_res slice_res 1]);
T = ([Mdc Pxyz0; 0 0 0 1]*D);

fclose(fid);
