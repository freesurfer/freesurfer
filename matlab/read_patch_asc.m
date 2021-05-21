function [surf] = read_patch_asc(fname)
% function [ surf] = read_asc(fname)
%
% Reads asc-file with patch-information


%
% read_asc.m
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

fp = fopen(fname, 'r');

% Dump first line
fgets(fp);

% Nr of vertices and faces
[a b] = fscanf(fp, '%d', 2);
nv = a(1) ;
nf = a(2) ;

fclose(fp);

% Read vertices and its indices
[ind vno x y z] = textread(fname, '%d vno=%d %f %f %f', nv, 'headerlines', 2, 'whitespace', '\t\b\n');

% Read faces and its indices
[find v1 v2 v3] = textread(fname, '%d %d %d %d', nf, 'headerlines', 2*nv+2, 'whitespace', '\t\b\n');

surf.faces = [v1 v2 v3];
surf.vertices = vno ;
surf.x = x;
surf.y = y;
surf.z = z;
