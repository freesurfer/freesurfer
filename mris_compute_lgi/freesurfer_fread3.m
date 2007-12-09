function [retval] = freesurfer_fread3(fid)

% freesurfer_fread3 - read a 3 byte integer out of a file
% 
% [retval] = freesurfer_fread3(fid)
%
% see also freesurfer_write3, freesurfer_read_surf, freesurfer_write_surf
% 

b1 = fread(fid, 1, 'uchar') ;
b2 = fread(fid, 1, 'uchar') ;
b3 = fread(fid, 1, 'uchar') ;
retval = bitshift(b1, 16) + bitshift(b2,8) + b3 ;