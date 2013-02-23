% You can place sources and destinations at any place of the mesh surface.
% Danil Kirsanov, 09/2007 

global geodesic_library;                
geodesic_library = 'geodesic_debug';      %"release" is faster and "debug" does additional checks
rand('state', 0);                         %comment this statement if you want to produce random mesh every time

N = 120;                                  %number of points in a mesh
[vertices,faces] = create_hedgehog_mesh(N, 0.1);   %create "noisy sphere" mesh; "vertices" contains 3D vertex coordinates; "faces" contains vertex id's for every triangle

[mesh, edge_to_vertex, edge_to_face] = geodesic_new_mesh(vertices,faces);         %initilize new mesh and receive edge info
disp(sprintf('mesh has %d edges', length(edge_to_vertex)));

source_points = {};
vertex_id = 10;         %put the first source at the vertex #10
source_points{1} = geodesic_create_surface_point('vertex',vertex_id,vertices(vertex_id,:));

face_id = 53;           %put the second source in the middle of the face #53
coords = mean(vertices(faces(face_id,:),:));        %compute 3D coordinates of this point
source_points{2} = geodesic_create_surface_point('face',face_id,coords);

edge_id = 85;           %put the third source on the edge #85
v1 = vertices(edge_to_vertex(edge_id,1),:);     %first vertex of the edge;
v2 = vertices(edge_to_vertex(edge_id,2),:);     %second vertex of the edge;
source_points{3} = geodesic_create_surface_point('edge',edge_id,0.2*v1+0.8*v2);

algorithm = geodesic_new_algorithm(mesh, 'exact');      %initialize new geodesic algorithm
geodesic_propagate(algorithm, source_points);   %propagation stage of the algorithm (the most time-consuming)

[source_ids, distances] = geodesic_distance_and_source(algorithm);     %for every vertex, figure out which source it belongs to

        %mesh is colored in three colors which show the 
hold off
trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3),source_ids, 'EdgeColor', 'k', 'FaceAlpha', 0.99, 'FaceColor', 'interp');       %plot the mesh
hold on;

colormap('spring');
c_map = colormap + 0.5;       % playing with colormap to improve visibility
colormap(c_map/max(c_map(:)));

destinations = {};
for vertex_id = 1:10:N;         %define destinations in a similar way
    destinations{end+1} = geodesic_create_surface_point('vertex',vertex_id,vertices(vertex_id,:));
end;

for face_id = 1:10:length(faces);           
    coeffs = rand(1,3);
    coeffs = coeffs/sum(coeffs);
    coords = coeffs*vertices(faces(face_id,:),:);    %put a point into a random place inside a face
    destinations{end+1} = geodesic_create_surface_point('face',face_id,coords);
end;

for edge_id = 1:10:length(edge_to_vertex);
    v1 = vertices(edge_to_vertex(edge_id,1),:);
    v2 = vertices(edge_to_vertex(edge_id,2),:);
    destinations{end+1} = geodesic_create_surface_point('edge',edge_id,0.4*v1+0.6*v2);
end

for i=1:length(destinations);
    path{i} = geodesic_trace_back(algorithm, destinations{i});     %find a shortest path from source to destination
    [x,y,z] = extract_coordinates_from_path(path{i});
    plot3(x,y,z);    %plot a sinlge path for this algorithm
    plot3(destinations{i}.x,destinations{i}.y,destinations{i}.z,'bo','MarkerSize',3, 'MarkerFaceColor', 'b');
end;

for i=1:length(source_points);
    plot3(source_points{i}.x,source_points{i}.y,source_points{i}.z,'ro','MarkerSize',5, 'MarkerFaceColor', 'r');
end;
daspect([1 1 1]);

geodesic_delete;