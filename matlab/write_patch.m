function write_patch(filename, patch)
% [curv] = write_curv(fname, curv, fnum)
%
% writes a patch into a binary file
%       filename - name of file to write to
%       patch    - patch structure
%

fid = fopen(filename, 'w', 'b');

% magic number
fwrite(fid, -1, 'int32');

% number of points in patch
fwrite(fid, patch.npts, 'int32');

% write each point
for i=1:patch.npts
    fwrite(fid, patch.ind(i) + 1, 'int32');
    fwrite(fid, patch.x(i), 'float');
    fwrite(fid, patch.y(i), 'float');
    fwrite(fid, patch.z(i), 'float');
end

fclose(fid);
