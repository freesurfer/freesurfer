function mris_write_m3z(filename_m3z, geom_orig, vol_orig, geom_dest, ...
                        vol_dest, vol_ind0, spacing, exp_k)
%MRIS_WRITE_M3Z
%
% mris_write_m3z(filename_m3z, orig_geom, vol_orig, dest_geom,
%                vol_dest, vol_ind0 [, spacing, exp_k])
%
%
% see also MRIS_READ_M3Z.
%
% see GCAMwrite() in gcamorph.c
%
% oliver hinds <ohinds@nmr.mgh.harvard.edu>, 2015/dec/29
%**************************************************************************%

    if ( nargin == 0 ), help(mfilename); return; end;

    if nargin < 7
        spacing = 1;
    end

    if nargin < 8
        exp_k = 0;
    end

    %==--------------------------------------------------------------------==%

    width = size(vol_orig, 1);
    height = size(vol_orig, 2);
    depth = size(vol_orig, 3);

    [fp, errstr] = fopen(filename_m3z, 'w', 'b');
    if ( fp == -1 ),
        error(errstr);
    end;

    % version
    fwrite(fp, [1], 'float');

    % dims, spacing
    fwrite(fp, [width height depth spacing], 'int');

    % exp_k
    fwrite(fp, [exp_k], 'float');

    %==--------------------------------------------------------------------==%

    % vectorized data write.
    nvox = width * height * depth;
    write_buf = repmat(char(0), 9 * 4 * nvox, 1);

    % reverse byte ordering, copy each volume into the buffer

    % orig
    float_buf = reshape(single(permute(vol_orig, [4 3 2 1])), 3 * nvox, 1);
    buf_offset = 9 * 4 * repmat(0:((nvox) - 1), 3 * 4, 1);
    inds = repmat([4:-1:1 8:-1:5 12:-1:9]', nvox, 1) + buf_offset(:);
    write_buf(inds) = typecast(float_buf, 'uint8');

    % dest
    float_buf = reshape(single(permute(vol_dest, [4 3 2 1])), 3 * nvox, 1);
    inds = inds + 12;
    write_buf(inds) = typecast(float_buf, 'uint8');

    % dest
    int_buf = reshape(int32(permute(vol_ind0, [4 3 2 1])), 3 * nvox, 1);
    inds = inds + 12;
    write_buf(inds) = typecast(int_buf, 'uint8');

    fwrite(fp, write_buf);

    % write vol geom info
    TAG_GCAMORPH_GEOM   = 10;
    fwrite(fp, [TAG_GCAMORPH_GEOM], 'int');

    write_geom(fp, geom_dest);
    write_geom(fp, geom_orig);

    fclose(fp);

    % hack.
    gzip(filename_m3z);
    movefile(strcat(filename_m3z, '.gz'), filename_m3z);

end

function write_geom(fp, geom)
    fwrite(fp, [1], 'int');
    fwrite(fp, geom.dim(2:4), 'int');
    fwrite(fp, geom.pixdim(2:4), 'float');

    for c=1:3
        for r=1:3
            fwrite(fp, [geom.vox2ras(r, c) / geom.pixdim(c + 1)], 'float');
        end
    end

    ras = geom.vox2ras * [geom.dim(2:4) ./ 2; 1];
    fwrite(fp, ras(1:3), 'float');

    % TODO write volume filenames
    fwrite(fp, zeros(512, 1, 'uint8'), 'uint8');
end
