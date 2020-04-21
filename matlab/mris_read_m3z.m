function [vol_orig, vol_dest, vol_ind0, spacing, exp_k] = mris_read_m3z(filename_m3z)
%MRIS_READ_M3Z
%
% varargout = mris_read_m3z(filename_m3z)
% function [vol_orig, vol_dest, vol_ind0, spacing, exp_k] = mris_read_m3z(filename_m3z)
%
%
% see also MRIS_SAVE_M3Z.

% see GCAMread() in gcamorph.c

% jonathan polimeni <jonp@nmr.mgh.harvard.edu>, 2012/feb/05
% oliver hinds <ohinds@nmr.mgh.harvard.edu>, 2015/dec/29
%**************************************************************************%

  if ( nargin == 0 ), help(mfilename); return; end;


  %==--------------------------------------------------------------------==%

  filename = gunzip(filename_m3z, tempname);
  fp = -1;

  % protect against dangling (potentially large) tempfiles
  try
    [fp, errstr] = fopen(filename{1}, 'r', 'b');
    if ( fp == -1 ),
      error(errstr);
    end;

    version = fread(fp, 1, 'float');

    if ( version ~= 1 ),
      error('currently only version 1.0 M3Z files are supported');
    end;

    width = fread(fp, 1, 'int');
    height = fread(fp, 1, 'int');
    depth = fread(fp, 1, 'int');
    spacing = fread(fp, 1, 'int');
    exp_k = fread(fp, 1, 'float');

    %==--------------------------------------------------------------------==%

    % vectorized data read. read all data into a buffer that can be
    % typecast into the appropriate data type all at once.

    % read all the data (9 numbers per voxel, each at 32 bits)
    buf = fread(fp, width * height * depth * 9 * 4, 'uchar=>uchar');

    % build a map of indices where the data resides.
    vox_offset = 9 * 4 * repmat(0:((width * height * depth) - 1), 3 * 4, 1);
    inds = repmat([4:-1:1 8:-1:5 12:-1:9]', width * height * depth, 1) + vox_offset(:);

    % extract the three interleaved volumes and permute the result to match the
    % original looped read. the double conversion isn't necessary, but is
    % included to maintain backward compatibility.
    vol_orig = reshape(typecast(buf(inds), 'single'), 3, depth, height, width);
    vol_orig = double(permute(vol_orig, [4 3 2 1]));

    inds = inds + 12;
    vol_dest = reshape(typecast(buf(inds), 'single'), 3, depth, height, width);
    vol_dest = double(permute(vol_dest, [4 3 2 1]));

    inds = inds + 12;
    vol_ind0 = reshape(typecast(buf(inds), 'int32'), 3, depth, height, width);
    vol_ind0 = double(permute(vol_ind0, [4 3 2 1]));

    fpos_dataend = ftell(fp);

    tag = fread(fp, 1, 'int');

    % jump to end of file
    fseek(fp, 0, 1);
    fpos_fileend = ftell(fp);

    databytes_remaining = fpos_fileend - fpos_dataend;

    disp(sprintf('==> [%s]: binary data remaining after displacement map: %10.2f MB', mfilename, databytes_remaining/2^20));

    % remaining data can represent TAG_GCAMORPH_LABELS, TAG_GCAMORPH_GEOM,
    % TAG_GCAMORPH_TYPE, or TAG_MGH_XFORM---see gcamorph.c

    TAG_GCAMORPH_GEOM   = 10;
    TAG_GCAMORPH_TYPE   = 11;
    TAG_GCAMORPH_LABELS = 12;
    TAG_MGH_XFORM       = 31;

    fclose(fp);

  catch
    delete(filename{1});
    rmdir(fileparts(filename{1}));

    display(lasterr);
    disp(err.stack);
  end

  delete(filename{1});
  rmdir(fileparts(filename{1}));
return;


  %************************************************************************%
  %%% $Source: /space/repo/1/dev/dev/matlab/mris_read_m3z.m,v $
  %%% Local Variables:
  %%% mode: Matlab
  %%% fill-column: 76
  %%% comment-column: 0
  %%% End:
