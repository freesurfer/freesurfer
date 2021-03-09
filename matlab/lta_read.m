function [M lta] = lta_read(fname)
% [M lta] = lta_read(fname)
% lta is a structure with elements from the LTA file
%  The lta structure is compatible with the LTA stucture in lta_write.m
%  

%
% lta_read.m
%
% Original Author: Bruce Fischl
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

if (strcmp(fname((length(fname)-2):length(fname)), 'xfm') | ...
	strcmp(fname((length(fname)-2):length(fname)), 'XFM'))
	M = xfm_read(fname) ;
	return ;
end


fid = fopen(fname) ;
if(fid < 0)
  error(sprintf('could not open file %s', fname));
end

% Go through the lines until get to to a line where the first
% string is the word 'type'
n=0;
while(1)
  n = n + 1;
  tline = fgetl(fid);
  if(tline == -1)
    fprintf('lta_read format error %s\n',fname);
    M = [];
    return;
  end
  tag = sscanf(tline,'%s',1);
  if(strcmp(tag,'type')) break; end
end
lta.type = sscanf(tline,'%*s %*s %d');

tline = fgetl(fid);   % nxforms
lta.nxforms = sscanf(tline,'%*s %*s %d');
tline = fgetl(fid);   % mean
lta.mean = sscanf(tline,'%*s %*s %f %f %f');
tline = fgetl(fid);   % sigma
lta.sigma = sscanf(tline,'%*s %*s %f');
tline = fgetl(fid);   % dimensions
lta.dims = sscanf(tline,'%d %d %d');

% Now read the matrix
M = zeros(4,4) ;
for row=1:4
  tline = fgetl(fid);   % one row of matrix
  tmp = sscanf(tline, '%f');
  M(row,:) = tmp';
end
if(nargout == 1) return; end

lta.xform = M;

tline = fgetl(fid);   % "src volume info"
tline = fgetl(fid);   % valid
lta.srcmri.valid = sscanf(tline,'%*s %*s %d');
tline = fgetl(fid);   % source filename
lta.srcfile = sscanf(tline,'%*s %*s %s');
tline = fgetl(fid);   % source volume dimensions
lta.srcmri.volsize = sscanf(tline,'%*s %*s %d %d %d');
tline = fgetl(fid);   % source voxel size
lta.srcmri.volres = sscanf(tline,'%*s %*s %f %f %f');
tline = fgetl(fid);   % xras
lta.srcmri.xras = sscanf(tline,'%*s %*s %f %f %f');
tline = fgetl(fid);   % yras
lta.srcmri.yras = sscanf(tline,'%*s %*s %f %f %f');
tline = fgetl(fid);   % zras
lta.srcmri.zras = sscanf(tline,'%*s %*s %f %f %f');
tline = fgetl(fid);   % cras
lta.srcmri.cras = sscanf(tline,'%*s %*s %f %f %f');
lta.srcmri.x_r = lta.srcmri.xras(1);
lta.srcmri.x_a = lta.srcmri.xras(2);
lta.srcmri.x_s = lta.srcmri.xras(3);
lta.srcmri.y_r = lta.srcmri.yras(1);
lta.srcmri.y_a = lta.srcmri.yras(2);
lta.srcmri.y_s = lta.srcmri.yras(3);
lta.srcmri.z_r = lta.srcmri.zras(1);
lta.srcmri.z_a = lta.srcmri.zras(2);
lta.srcmri.z_s = lta.srcmri.zras(3);
lta.srcmri.c_r = lta.srcmri.cras(1);
lta.srcmri.c_a = lta.srcmri.cras(2);
lta.srcmri.c_s = lta.srcmri.cras(3);

tline = fgetl(fid);   % "dst volume info"
tline = fgetl(fid);   % valid
lta.dstmri.valid = sscanf(tline,'%*s %*s %d');
tline = fgetl(fid);   % destination filename
lta.dstfile = sscanf(tline,'%*s %*s %s');
tline = fgetl(fid);   % destination volume dimensions
lta.dstmri.volsize = sscanf(tline,'%*s %*s %d %d %d');
tline = fgetl(fid);   % destination voxel size
lta.dstmri.volres = sscanf(tline,'%*s %*s %f %f %f');
tline = fgetl(fid);   % xras
lta.dstmri.xras = sscanf(tline,'%*s %*s %f %f %f');
tline = fgetl(fid);   % yras
lta.dstmri.yras = sscanf(tline,'%*s %*s %f %f %f');
tline = fgetl(fid);   % zras
lta.dstmri.zras = sscanf(tline,'%*s %*s %f %f %f');
tline = fgetl(fid);   % cras
lta.dstmri.cras = sscanf(tline,'%*s %*s %f %f %f');
lta.dstmri.x_r = lta.dstmri.xras(1);
lta.dstmri.x_a = lta.dstmri.xras(2);
lta.dstmri.x_s = lta.dstmri.xras(3);
lta.dstmri.y_r = lta.dstmri.yras(1);
lta.dstmri.y_a = lta.dstmri.yras(2);
lta.dstmri.y_s = lta.dstmri.yras(3);
lta.dstmri.z_r = lta.dstmri.zras(1);
lta.dstmri.z_a = lta.dstmri.zras(2);
lta.dstmri.z_s = lta.dstmri.zras(3);
lta.dstmri.c_r = lta.dstmri.cras(1);
lta.dstmri.c_a = lta.dstmri.cras(2);
lta.dstmri.c_s = lta.dstmri.cras(3);

lta.subject = '';
lta.fscale = .15;
tline = fgetl(fid); % subject, if there
while(tline ~= -1)
  tmp = sscanf(tline,'%s %*s');
  if(strcmp(tmp,'subject'))
    lta.subject = sscanf(tline,'%*s %s');
  end
  if(strcmp(tmp,'fscale'))
    lta.fscale = sscanf(tline,'%*s %f');
  end
  tline = fgetl(fid); % subject, if there
end

% Compute the vox2ras matrix 
Mdc = [lta.srcmri.xras lta.srcmri.yras lta.srcmri.zras];
Nvox2 = lta.srcmri.volsize/2;
D = diag(lta.srcmri.volres);
P0 = lta.srcmri.cras - Mdc*D*Nvox2;
lta.srcmri.vox2ras0 = [Mdc*D P0; 0 0 0 1];

Mdc = [lta.dstmri.xras lta.dstmri.yras lta.dstmri.zras];
Nvox2 = lta.dstmri.volsize/2;
D = diag(lta.dstmri.volres);
P0 = lta.dstmri.cras - Mdc*D*Nvox2;
lta.dstmri.vox2ras0 = [Mdc*D P0; 0 0 0 1];


fclose(fid) ;

