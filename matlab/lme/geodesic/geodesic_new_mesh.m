function [mesh, edge_to_vertex, edge_to_face] = geodesic_new_mesh(points, tri)

global geodesic_library
if ~libisloaded(geodesic_library)
    hfile = 'geodesic_matlab_api.h';
    loadlibrary([geodesic_library '.so'], hfile);
end

dim = find(size(points) == 3);
if dim == 1
    p = points(:);
else 
    p = points';
    p = p(:);
end;

dim = find(size(tri) == 3);
if dim == 1
    t = tri(:) - 1;
else 
    t = tri';
    t = t(:) - 1;
end;

tmp = libpointer('doublePtrPtr');
[id, tmp1, tmp2, num_edges, edges] = calllib(geodesic_library, 'new_mesh', length(p)/3, p, length(t)/3, t, 1, tmp);
setdatatype(edges, 'doublePtr', 4, num_edges);

mesh.id = id;
mesh.object_type = 'mesh';
edge_to_vertex = (edges.Value(1:2,:) + 1)';
edge_to_face = (edges.Value(3:4,:) + 1)';