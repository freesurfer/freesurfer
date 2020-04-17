function Surf = fs_read_surf(fname)
% Surf = fs_read_surf(fname)
% 
% Reads coordinates and triangles of a Freesurfer's surf (Depends on 
% Freesurfer's read_surf function).
%
% Input
% fname: A Freesurfer's surface file.
%
% Output
% Surf: Surface which represent a coordinate system where the analysis is 
% to be done. It is a structure with Surf.tri = t x 3 matrix of triangle 
% indices, 1-based, t=#triangles and Surf.coord = 3 x nv matrix of 
% coordinates, nv=#vertices.
%
% Original Author: Jorge Luis Bernal Rusiel 
%
[coord,faces] = read_surf(fname);
Surf.tri = uint32(faces + 1);
Surf.coord = coord';
