% Simple application: drawing approximate piece-wise linear isodistance curves (isolines).
% In theory, isolines on a triangular mesh are piecewise second-order curves.
% Here, we draw piecewise linear approximation of these curves. The distance function 
% of the mesh is sampled at the subdivision vertices within each triangle and 
% represented as a linear function within each subdivision triangle.
% This particular method works well only for 'exact' geodesic method,
% because for the approximate methods the distance field is discontinuous
% on the edges of the mesh.
% Danil Kirsanov, 09/2007 

subdivision_level = 5;     %number of additional vertices per edge; each mesh triangle will have roughly speaking subdivision_level^2 number of subdivision triangles
number_of_isolines = 10;                     

global geodesic_library;                
geodesic_library = 'geodesic_debug';      %"release" is faster and "debug" does additional checks
rand('state', 0);                         %comment this statement if you want to produce random mesh every time

[vertices,faces] = create_flat_triangular_mesh(0.2, 0); 
N = length(vertices);
% N = 300;                                  %number of points in a mesh
% [vertices,faces] = create_hedgehog_mesh(N, 0.1);   %create "noisy sphere" mesh; "vertices" contains 3D vertex coordinates; "faces" contains vertex id's for every triangle

hold off;
trisurf(faces,vertices(:,1),vertices(:,2),vertices(:,3), 'FaceColor', 'w', 'EdgeColor', 'k', 'FaceAlpha', 0.99);       %plot the mesh
hold on;

mesh = geodesic_new_mesh(vertices,faces);         %initilize new mesh
algorithm = geodesic_new_algorithm(mesh, 'exact');      %initialize new geodesic algorithm

source_points = {};
for vertex_id = [26,N-5, floor(N/2)-3, floor(N/2)+3];          %create source points 
    source_points{end+1} = geodesic_create_surface_point('vertex',vertex_id, vertices(vertex_id,:));
end;

geodesic_propagate(algorithm, source_points); 

[sub_weights,sub_tri] = create_subdivision_pattern(subdivision_level); %coordinates of the subdivision triangles

[source_id, distances] = geodesic_distance_and_source(algorithm);     %find distances for all vertices of the mesh
max_distance = max(distances);

isolines_step = max_distance/number_of_isolines;
isolines_distances = isolines_step/2 : isolines_step : max_distance;

[known_distances, tmp] = find(sub_weights > 0.999);                %we already know the distances to the actual vertices of the mesh
unknown_distances = setdiff([1:size(sub_weights,1)], known_distances);
for face_id = 1:size(faces,1)      
    vertex_ids = faces(face_id,:);
    sub_vertices = sub_weights*vertices(vertex_ids,:);         %find coordinates of subdivision vertices
    sub_distances = zeros(size(sub_vertices,1),1);
    
    sub_distances(known_distances) = distances(vertex_ids);     %no need to reevaluate distance finction for the vertices of the face
    %for sub_index = 1:size(sub_vertices,1)                  %find distances at all subdivision vertices
    for sub_index = unknown_distances          %find distances at all internal subdivision vertices
        q = geodesic_create_surface_point('face',face_id,sub_vertices(sub_index,:));
        [source_id, d] = geodesic_distance_and_source(algorithm, q);
        sub_distances(sub_index) = d;
    end;
    
    selected_distances = find(isolines_distances >= min(sub_distances) & ...
                              isolines_distances <= max(sub_distances));
    selected_distances = isolines_distances(selected_distances);
                          
    for tri_index = 1:size(sub_tri,1)               %intersection of the equidistance curve with every subdivision triangle is approximated as a straight line
        for isodistance = selected_distances
            p = [];
            for i = 1:3
                index1 = sub_tri(tri_index, i);
                index2 = sub_tri(tri_index, mod(i,3) + 1);
                d1 = sub_distances(index1);
                d2 = sub_distances(index2);
                
                if d1 == isodistance        %special case when isodistance paths exactly through the vertex
                    p(end+1,:) = sub_vertices(index1,:);
                end
                
                if isodistance > min(d1,d2) & isodistance < max(d1,d2)   %draw equidistant curve
                    a = (isodistance - d2)/(d1 - d2);       % I believe that division by zero will never happen here
                    p(end+1,:) = a*sub_vertices(index1,:) + (1-a)*sub_vertices(index2,:);
                end;
            end
            if size(p,1) == 2
                 plot3(p(:,1), p(:,2),p(:,3),'b', 'LineWidth', 2);
            end
        end
    end;
end;

for i=1:length(source_points);
    h = plot3(source_points{i}.x,source_points{i}.y,source_points{i}.z,'ro','MarkerSize',5, 'MarkerFaceColor', 'r');
end;

campos([0, 0, 5])
legend(h, 'sources');

geodesic_delete;