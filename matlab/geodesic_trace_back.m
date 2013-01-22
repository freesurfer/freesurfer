function path = geodesic_trace_back(algorithm, destination)

global geodesic_library;

tmp{1} = destination;
d = geodesic_convert_surface_points(tmp);

tmp = libpointer('doublePtrPtr');
[path_length, tmp, path_double] = calllib(geodesic_library, 'trace_back', algorithm.id, d, tmp);

setdatatype(path_double, 'doublePtr', path_length*5);
path = geodesic_convert_surface_points(path_double.Value);