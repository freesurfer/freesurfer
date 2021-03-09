function m = fast_ldbhdr(bhdrfile);
%
% m = fast_ldbhdr(bhdrfile)
%
% Parses the bhdr header file for bshorts/bfloats and returns the
% result as an mri structure. See fast_mri_struct. The bhdr file
% consists of keyword/value pairs.
%
% bhdrfile is the name of the bhdr file. If this cannot be opened,
% then bhdrfile.bhdr is tried (ie, it treats bhdrfile as the stem).
% If this cannot be opened, then it retuns null.
%
% top_left_?:     RAS at Top-Left  Edge = T*[0 0 0 1]'
% top_right_?:    RAS at Top-Right Edge = T*[ncols 0 0 1]'
% bottom_right_?: RAS at Bot-Right Edge = T*[ncols nrows 0 1]'
%
% See also fast_svbhdr and fast_mri_struct.
%
%


%
% fast_ldbhdr.m
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

m = [];

if(nargin ~= 1)
  msg = 'm = fast_ldbhdr(bhdrfile)';
  qoe(msg); error(msg);
end

fid = fopen(bhdrfile,'r');
if(fid == -1)
  bhdrfilebhdr = sprintf('%s.bhdr',bhdrfile);
  fid = fopen(bhdrfilebhdr,'r');
  if(fid == -1) return; end    
end

m = fast_mri_struct;

ncols = [];

%------- Loop through all the lines -----------%
nthrow = 1;
while(1)

  % scroll through any blank lines or comments %
  while(1)
    tline = fgetl(fid);
    if(~isempty(tline) & tline(1) ~= '#') break; end
  end
  if(tline(1) == -1) break; end

  % Read the keyword and value %  
  key = sscanf(tline,'%s',1);
  val = sscanf(tline,'%*s %f',1);
  %fprintf('key = %s, val = %g\n',key,val);
  
  % determine where to put it %
  switch(key)
   case {'cols:'},          ncols = val;
   case {'rows:'},          nrows = val;
   case {'nslices:'},       nslices = val;
   case {'n_time_points:'}, m.nframes = val;
   case {'slice_thick:'},   thickness = val;
   case {'top_left_r:'},    tlr = val;
   case {'top_left_a:'},    tla = val;
   case {'top_left_s:'},    tls = val;
   case {'top_right_r:'},   trr = val;
   case {'top_right_a:'},   tra = val;
   case {'top_right_s:'},   trs = val;
   case {'bottom_right_r:'},brr = val;
   case {'bottom_right_a:'},bra = val;
   case {'bottom_right_s:'},brs = val;
   case {'normal_r:'},      nr = val;
   case {'normal_a:'},      na = val;
   case {'normal_s:'},      ns = val;
   case {'image_te:'},      m.te = val;
   case {'image_tr:'},      m.tr = 1000*val; % msec
   case {'image_ti:'},      m.ti = val;
   case {'flip_angle:'},    m.flip_angle = val;
  end
  
end % while (1)
fclose(fid);
%---------------------------------------------%

if(isempty(ncols))
  fprintf('ERROR: fast_ldbhdr: reading %s\n',bhdrfile);
  m = [];
  return;
end
m.voldim = [ncols nrows nslices]'; %'

TL = [tlr tla tls]';
TR = [trr tra trs]';
BR = [brr bra brs]';
m.sdc = [nr na ns]'; 

% Compute the Column and Row FOV %
% Outer edge of voxel to outer edge of voxel
col_fov = sqrt(sum((TR-TL).^2));
row_fov = sqrt(sum((TR-BR).^2));

if(col_fov == 0 | row_fov == 0)
  % Probably a surface %
  return;
end

% Compute the Column and Row Resolution %
col_res = col_fov/ncols;
row_res = row_fov/nrows;

% Compute the Column and Row Normal Vectors %
m.cdc = (TR-TL)/col_fov;
m.rdc = (BR-TR)/row_fov;

m.volres = [col_res row_res thickness]';%'
m.P0 = TL; % Assumes TL is center, not edge
Mdc = [m.cdc m.rdc m.sdc];
D = diag(m.volres);
m.T = [Mdc*D m.P0; 0 0 0 1];

% Fix P0 for when TL, TR, and BR are the edges
%tmp = m.T*[0.5 0.5 0.5 1]';
%m.P0 = tmp(1:3);
%m.T(:,4) = tmp;

m.c = m.T*[m.voldim/2; 1];

return;
