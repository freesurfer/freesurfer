function err = fast_svbhdr(m, bhdrfile)
% err = fast_svbhdr(m, bhdrfile)
%
% m is a fast_mri_struct.
% It is assumed that zero-based indexing is used.
%
% See also fast_ldbhdr and fast_mri_struct.
%
% $Id: fast_svbhdr.m,v 1.2 2003/08/02 00:56:14 greve Exp $

err = 1;

if(nargin ~= 2)
  fprintf('err = fast_svbhdr(fsfmri, bhdrfile)\n');
  return;
end

fid = fopen(bhdrfile,'w');
if(fid == -1)
  fprintf('ERROR: could not open %s\n',bhdrfile);
  return;
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
fprintf(fid,'      image_tr: %f\n',m.tr);
fprintf(fid,'      image_ti: %f\n',m.ti);
fprintf(fid,'    flip_angle: %f\n',pi*m.flip_angle/180);

fclose(fid);

err = 0;

return;
