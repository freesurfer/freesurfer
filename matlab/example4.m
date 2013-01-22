% Advanced tuning: to save time, you can limit the propagation.
% Danil Kirsanov, 09/2007 

global geodesic_library;                
geodesic_library = 'geodesic_debug';      %"release" is faster and "debug" does additional checks
rand('state', 0);                         %comment this statement if you want to produce random mesh every time

[vertices,faces] = create_flat_triangular_mesh(0.2, 0); 
N = length(vertices);

mesh = geodesic_new_mesh(vertices,faces);         %initilize new mesh
algorithm = geodesic_new_algorithm(mesh, 'exact');      %initialize new geodesic algorithm

source_points = {};
for vertex_id = [26,N-5];          %create two source points 
    source_points{end+1} = geodesic_create_surface_point('vertex',vertex_id, vertices(vertex_id,:));
end;

stop_points = {};
for vertex_id = [29,N-14];          %create stop points
    stop_points{end+1} = geodesic_create_surface_point('vertex',vertex_id, vertices(vertex_id,:));
end;

stop_distance = 0.1;            %define stop distance

% propagation stops when covers stop_points and exceeds stop_distance
% it saves time and memory, but you have to be careful working with results
geodesic_propagate(algorithm, source_points, stop_points, stop_distance); 

distances = zeros(N,1);     %for every vertex, figure out which source it belongs to
paths = {};
for i=1:N   
    q = geodesic_create_surface_point('vertex',i,vertices(i,:));   %create destination at vertex #i
    [source_id, d] = geodesic_distance_and_source(algorithm, q);    
    distances(i) = d;
    path{i} = {};
    if(d < 1e10)        %if distance is found, we can trace a path
        path{i} = geodesic_trace_back(algorithm, q);     %find a shortest path from source to destination
    end
end;

undefined = find(distances > 1e10);     %find vertices for which the distances are undefined
defined = find(distances < 1e10);     
distances(undefined) = max(distances(defined))*1.1;   %plot undefined distances with white
colormap('winter');
c_map = colormap;       
c_map(end,:) = 1;
colormap(c_map);

hold off;
trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3),distances, 'FaceColor', 'interp', 'EdgeColor', 'k');       %plot the mesh
hold on;
daspect([1 1 1]);

for i=1:length(path)                   %plot paths to all vertices with finite distances                                         
    [x,y,z] = extract_coordinates_from_path(path{i});
    plot3(x,y,z,'r');    
end;

for i=1:length(stop_points);        %show stop vertices
    h = plot3(stop_points{i}.x,stop_points{i}.y,stop_points{i}.z,'ro','MarkerSize',5, 'MarkerFaceColor', 'r');
end;
legend(h, 'stop vertices');

geodesic_delete;
