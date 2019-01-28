function [vtxs,nvtxs] = fs_read_label(labelf)
% [vtxs,nvtxs] = fs_read_label(labelf)
% 
% Reads the indices of the vertices of a Freesurfer's label. 
%
% Input
% labelf: A Freesurfer's label file.
%
% Output
% vtxs: Indices of the vertices (1-based). 
% nvtxs: The number of vertices in the label.
%
% $Revision: 1.2 $  $Date: 2015/01/06 17:14:50 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: mreuter $
%    $Date: 2015/01/06 17:14:50 $
%    $Revision: 1.2 $
%
fid = fopen(labelf);
tline = fgetl(fid);
nvtxs = fscanf(fid,'%d',1);
vtxs = uint32(fscanf(fid,'%d %*g %*g %*g %*g',[1 nvtxs])); 
vtxs = vtxs + 1;
fclose(fid); 

