function m = fast_mri_struct
% create a default structure for an mri. Note: all values
% assume zero-based indexing.
%
% $Id: fast_mri_struct.m,v 1.1 2003/08/02 00:58:37 greve Exp $

m.te         = 0; % msec
m.tr         = 0; % sec
m.ti         = 0; % msec
m.flip_angle = 0; % degrees

m.voldim = [0 0 0]; % ncols nrows nslices
m.nframes = [0];    % number of frames, time points, planes, etc
m.T = zeros(4);     % zero-based vox2ras

% These can be derived from the above
m.volres = [0 0 0]; % dcol, drow, dslice
m.cdc = [];    % column direction cos
m.rdc = [];    % row direction cos
m.sdc = [];    % slice direction cos
m.P0  = [];    % position at center of voxel 0,0,0
m.c   = [];    % The "center" of the volume as def by FreeSurfer

return;







