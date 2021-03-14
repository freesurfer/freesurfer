function [S, v, f, ORIGIN] = read_vf(fname, patch)
% function [S, v, f, ORIGIN] = read_vf(fname, patch)
% This function reads a matfile with vertices/faces-information in it.
% if [fname '.mat'] cannot be found, a '.geo'- or 'asc'-file is searched
% converted and saved as a matfile.
% S(1) = nr of surfaces
% S(2) = nr of vertices
% S(3) = nr of faces
% v: n x 3 matrix of vertices
% f: m x 4 matrix of faces
% ORIGIN: coordinates (2D or 3D) of ORIGIN in mm space, where the coordinates
% [0 0 0] in vertices space are mapped to.

% fname: filename of mat-file to be opened
% patch: flag: if 1, try to get patch-information, i.e. indices of vertices and faces
% in generating original surface, if 0, we don't need this information


%
% read_vf.m
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

fn = spm_str_manip(fname, 'r');

source = 0;

if ~exist([fn '.mat'])
   % Try to open asc-file
   if patch == 0 | ~exist([fn '.asc'], 'file')
      % Try to open geo-file
      if patch == 1 
         % no mat- or asc-file
         error(sprintf('Cannot find patch-infomation, file %s not found', [fn '.asc']))
      elseif ~exist([fn '.geo'], 'file')
         error(sprintf('Cannot find %s or convert a geo- or asc-file to it.', [fn '.mat']))
      else
         [S, v, f] = read_moviebyu([fn '.geo']);
         source = 2;
      end
   else
      [S, v, fi] = read_asc([fn '.asc']);
      source = 1;
   end
else
   load([fn '.mat']);
   if patch == 1 & exist('fi') ~= 1
      % used asked for patch information, but mat-file didn't contain this information
      % so open asc-file (can happen, if asc-file wasn't present when mat-file
      % was generated
      if exist ([fn '.asc'], 'file')
         [S, v, fi] = read_asc([fn '.asc']);
         save([fn '.mat'], 'fi', 'v', '-append');
         source = 1;
      else
         error(sprintf('Cannot find %s', [fn '.asc']));
      end
   elseif patch == 0 & exist('f') ~= 1
      % no patch information needed, try to get plain f
      if exist ([fn '.geo'], 'file')
         [S, v, f] = read_moviebyu([fn '.geo']);
         save([fn '.mat'], 'f', '-append');
         source = 2;
      else
         error(sprintf('Cannot find %s', [fn '.geo']));
      end
   end
   
   if patch == 0 & size(v, 2) ~= 3
      v = v(:, 2:4);
   end
end

if source == 1 | source == 2
   % Assume that ORIGIN is [128 128 128]
   ORIGIN = [128 128 128];
end

if source == 1
   save([fn '.mat'], 'S', 'v', 'fi', 'ORIGIN');
elseif source == 2
   save([fn '.mat'], 'S', 'v', 'f', 'ORIGIN'); 
end

if patch == 1
   f = fi;
end   





































































































