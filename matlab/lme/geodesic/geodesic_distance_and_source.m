%finds best source and distance to the best source
% if distance is negative, the best source cannot be found (for example, because the propagation was stopped before it reached this point)
% Danil Kirsanov, 09/2007 

function [source_id, distance] = geodesic_best_source(algorithm, destination)

global geodesic_library;

if nargin == 2
    d = geodesic_convert_surface_points({destination});

    tmp = 1;
    [source_id, tmp1, distance] = calllib(geodesic_library, 'distance_and_source', algorithm.id, d, tmp);
    source_id = source_id + 1;
else                                    %return distances and sources for all vertices
    tmp = libpointer('doublePtrPtr');
%     tmp = libpointer('int32PtrPtr');
     tmp1 = libpointer('longPtrPtr');
%  tmp1 = libpointer('doublePtr');
    
    [num_vertices, d, s] = calllib(geodesic_library, 'distance_and_source_for_all_vertices', algorithm.id, tmp, tmp1);
    
    setdatatype(d, 'doublePtr', num_vertices);
    distance = d.Value;
    setdatatype(s, 'longPtr', num_vertices);
    source_id = s.Value + 1;
end