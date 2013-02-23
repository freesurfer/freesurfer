function geodesic_propagate(algorithm, source_points, stop_points, max_distance)

global geodesic_library;

if nargin < 4
    max_distance = 1e100;
end

if nargin < 3
    stop_points = [];
end

if ~libisloaded(geodesic_library)
    disp('error: geodesic library is not loaded');
    return;
end

sources = geodesic_convert_surface_points(source_points);
stops = geodesic_convert_surface_points(stop_points);

calllib(geodesic_library, 'propagate', algorithm.id, ...
    sources, length(sources)/5, stops, length(stops)/5, max_distance);
