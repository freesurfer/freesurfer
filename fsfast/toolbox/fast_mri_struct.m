function m = fast_mri_struct(orientation)
% m = fast_mri_struct(<orientation>)
%
% Create a default structure for an mri. Note: all values
% assume zero-based indexing. Without any inputs, none of 
% the values are set.
%
% If The direction cosines can be set with orientation:
%   LRAPSI
%   LRAPIS
%   LRPASI
%   LRPAIS
%
% $Id: fast_mri_struct.m,v 1.2 2004/10/06 20:46:30 greve Exp $

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

if(~exist('orientation','var')) return; end

switch(orientation)
 case 'LRAPSI'
  m.cdc = [-1  0  0]';
  m.rdc = [ 0 -1  0]';
  m.sdc = [ 0  0 -1]';
 case 'LRAPIS'
  m.cdc = [-1  0  0]';
  m.rdc = [ 0 -1  0]';
  m.sdc = [ 0  0 +1]';
 case 'LRPASI'
  m.cdc = [-1  0  0]';
  m.rdc = [ 0 +1  0]';
  m.sdc = [ 0  0 -1]';
 case 'LRPAIS'
  m.cdc = [-1  0  0]';
  m.rdc = [ 0 +1  0]';
  m.sdc = [ 0  0 +1]';
 otherwise
  fprintf('ERROR: unrecognized orientation %s\n',orientation);
end


return;







