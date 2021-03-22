function err = fast_svbhdr(m, bhdrfile, isstem)
% err = fast_svbhdr(m, bhdrfile, <isstem>)
%
% m is a fast_mri_struct.
% bhdrfile is the name of the bhdrfile or a bstem, depending
%  upon the presence and value of isstem.
% isstem - if present and equal to 1 then bhdrfile treated as
%  a stem and the bhdr file name is bhdrfile.bhdr
%
% It is assumed that zero-based indexing is used.
%
% See also fast_ldbhdr and fast_mri_struct.
%
%


%
% fast_svbhdr.m
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

err = 1;

if(nargin ~= 2 & nargin ~= 3)
  fprintf('err = fast_svbhdr(fsfmri, bhdrfile, <isstem>)\n');
  return;
end

if(exist('isstem') ~= 1) isstem = []; end
if(isempty(isstem)) isstem = 0; end

if(isstem)
  bhdrfile = sprintf('%s.bhdr',bhdrfile);
end

fid = fopen(bhdrfile,'w');
if(fid == -1)
  fprintf('ERROR: could not open %s\n',bhdrfile);
  return;
end
  
if(mar(m.T)==0 & ~isempty(m.cdc))
  % Does not look like vox2ras is set, so create from dircos
  Mdc = [reshape1d(m.cdc) reshape1d(m.rdc) reshape1d(m.sdc)];
  D = diag(m.volres);
  if(isempty(m.P0)) P0 = [0 0 0]'; end
  m.T = [Mdc*D reshape1d(P0); 0 0 0 1];
end


TL = m.T*[0 0 0 1]'; % Center of first vox
TR = m.T*[m.voldim(1) 0 0 1]'; % Top-right edge + 0.5 vox
BR = m.T*[m.voldim(1) m.voldim(2) 0 1]'; % Bot right + 0.5 vox

% Valid for when TL, TR, and BR are the edges
%TL = m.T*[-0.5 -0.5 -0.5 1]'; % Top-left edge
%TR = m.T*[m.voldim(1)-0.5 -0.5 -0.5 1]'; % Top-right edge
%BR = m.T*[m.voldim(1)-0.5 m.voldim(2)-0.5 -0.5 1]'; % Bot right

SliceNorm = m.T(:,3);
mag = sqrt(sum(SliceNorm.^2));
if(mag>0) SliceNorm = SliceNorm/mag; end

fprintf(fid,'          cols: %d\n',m.voldim(1));
fprintf(fid,'          rows: %d\n',m.voldim(2));
fprintf(fid,'       nslices: %d\n',m.voldim(3));
fprintf(fid,' n_time_points: %d\n',m.nframes);
fprintf(fid,'   slice_thick: %d\n',m.volres(3));
fprintf(fid,'    top_left_r: %f\n',TL(1));
fprintf(fid,'    top_left_a: %f\n',TL(2));
fprintf(fid,'    top_left_s: %f\n',TL(3));
fprintf(fid,'   top_right_r: %f\n',TR(1));
fprintf(fid,'   top_right_a: %f\n',TR(2));
fprintf(fid,'   top_right_s: %f\n',TR(3));
fprintf(fid,'bottom_right_r: %f\n',BR(1));
fprintf(fid,'bottom_right_a: %f\n',BR(2));
fprintf(fid,'bottom_right_s: %f\n',BR(3));
fprintf(fid,'      normal_r: %f\n',SliceNorm(1));
fprintf(fid,'      normal_a: %f\n',SliceNorm(2));
fprintf(fid,'      normal_s: %f\n',SliceNorm(3));
fprintf(fid,'      image_te: %f\n',m.te);
fprintf(fid,'      image_tr: %f\n',m.tr/1000); %msec->sec
fprintf(fid,'      image_ti: %f\n',m.ti);
fprintf(fid,'    flip_angle: %f\n',m.flip_angle);

fclose(fid);

err = 0;

return;
