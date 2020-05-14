function [vertices, faces] = freesurfer_read_surf(fname)

% freesurfer_read_surf - FreeSurfer I/O function to read a surface file
% 
% [vertices, faces] = freesurfer_read_surf(fname)
% 
% Reads the vertex coordinates (mm) and face lists from a surface file.
% 
% Surface files are stored as either triangulations or quadrangulations.
% That is, for a triangulation, each face is defined by 3 vertices.  For a
% quadrangulation, each face is defined by 4 vertices.  The rows of 'faces'
% contain indices into the rows of 'vertices', the latter holds the XYZ
% coordinates of each vertex.
%
% The freesurfer faces index the vertices in counter-clockwise order (when
% viewed from the outside of the surface).  This is consistent with a
% right-hand rule.  If we have vertices
%
% C           B
%
%
%       A
%
% Then we can calculate an edge vector from A to B (ie, AB = B - A) and
% another edge vector from A to C (ie, AC = C - A).  If you form a "gun"
% with your thumb and forefinger of the right hand, then align your thumb
% with the AB vector and your forefinger with the AC vector, your palm is
% facing out of the screen and extending your middle finger in the
% orthogonal direction to the plane of the screen will give the outward
% surface normal of the triangle ABC.  (If you lookup "triangle" on
% Wolfram's mathworld, you can see that AB is referred to as c and AC is
% referred to as b.)
%
% However, if this surface is read into matlab, it will give INWARD surface
% normals in the matlab patch command.  For some reason, matlab is not
% following the right hand rule.  To get OUTWARD normals with the matlab
% patch command, use faces(:,[1 3 2]) (see below).
%
% The vertex coordinates are in mm.  The FreeSurfer coordinate
% system for surfaces is quite simple, but relating to their MRI
% cor-??? files is too confusing to explain here; see the FreeSurfer
% homepage or google the documentation by Graham Wideman.  For the
% surfaces, at least, the origin is somewhere in the center of the
% head, and the vertex XYZ coordinates are oriented such that +X is
% right, +Y is anterior and +Z is superior (this is the
% FreeSurfer RAS coordinate system).
%
% Note that reading the faces of a quad file can take a long
% time due to their compact storage format.  In this case, the return of
% vertices can be faster if the face output variable is not specified; in
% this case, the faces are not read.
% 
% Try this to visualize the surface:
% Hp = patch('vertices',vertices,'faces',faces(:,[1 3 2]),...
%       'facecolor',[.5 .5 .5],'edgecolor','none')
% camlight('headlight','infinite')
% vertnormals = get(Hp,'vertexnormals');
%
% See also freesurfer_write_surf, freesurfer_read_curv,
%          freesurfer_read_wfile
%


% Copyright (C) 2000  Darren L. Weber

% History:  08/2000, Darren.Weber_at_radiology.ucsf.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% There used to be a line here where ver was set to the version of
% the file. Without it, "ver" reverts to a matlab defined structure
% and generates an error.
%fprintf('FREESURFER_READ_SURF [v %s]\n',ver(11:15));

if(nargin < 1)
    help freesurfer_read_surf;
    return;
end

%QUAD_FILE_MAGIC_NUMBER =  (-1 & 0x00ffffff) ;
%NEW_QUAD_FILE_MAGIC_NUMBER =  (-3 & 0x00ffffff) ;

TRIANGLE_FILE_MAGIC_NUMBER  =  16777214;
QUAD_FILE_MAGIC_NUMBER      =  16777215;


% open it as a big-endian file
fid = fopen(fname, 'rb', 'b');
if (fid < 0),
    str = sprintf('could not open surface file %s.', fname);
    error(str);
end

fprintf('...reading surface file: %s\n', fname);
tic;

magic = freesurfer_fread3(fid);

if (magic == QUAD_FILE_MAGIC_NUMBER),
    Nvertices = freesurfer_fread3(fid);
    Nfaces = freesurfer_fread3(fid);
    fprintf('...reading %d quad file vertices\n',Nvertices);
    vertices = fread(fid, Nvertices*3, 'int16') ./ 100 ; 
    if (nargout > 1),
        fprintf('...reading %d quad file faces (please wait)\n',Nfaces);
        faces = zeros(Nfaces,4);
        for iface = 1:Nfaces,
            for n=1:4,
                faces(iface,n) = freesurfer_fread3(fid) ;
            end
            if(~rem(iface, 10000)), fprintf(' %7.0f',iface); end
            if(~rem(iface,100000)), fprintf('\n'); end
        end
    end
elseif (magic == TRIANGLE_FILE_MAGIC_NUMBER),
    fprintf('...reading triangle file\n');
    tline = fgets(fid); % read creation date text line
    tline = fgets(fid); % read info text line
    
    Nvertices = fread(fid, 1, 'int32'); % number of vertices
    Nfaces = fread(fid, 1, 'int32'); % number of faces
    
    % vertices are read in column format and reshaped below
    vertices = fread(fid, Nvertices*3, 'float32');
    
    % faces are read in column format and reshaped
    faces = fread(fid, Nfaces*3, 'int32');
    faces = reshape(faces, 3, Nfaces)';
else
    str = sprintf('unknown magic number in surface file %s.', fname);
    error(str);
end

vertices = reshape(vertices, 3, Nvertices)';
fclose(fid);

fprintf('...adding 1 to face indices for matlab compatibility.\n');
faces = faces + 1;

t=toc; fprintf('...done (%6.2f sec)\n\n',t);

return
