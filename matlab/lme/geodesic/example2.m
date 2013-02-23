% This example shows how to work with multiple meshes and algorithms simultaneously.
% Danil Kirsanov, 09/2007 

global geodesic_library;                
geodesic_library = 'libgeodesic_debug';      %"release" is faster and "debug" does additional checks

rand('state', 0);                         %comment this statement if you want to produce random mesh every time
analysis = {};
                    
% for simplicity, we are going to keep all data related to
% the first and second meshes in analysis{1} and analysis{2} respectively

analysis{1}.N = 300;         %generate first mesh
[analysis{1}.vertices,analysis{1}.faces] = create_hedgehog_mesh(analysis{1}.N, 0.1, 0);   %create fist 

[analysis{2}.vertices,analysis{2}.faces] = create_flat_triangular_mesh(0.2, 0); 
analysis{2}.N = length(analysis{2}.vertices);  %another simple mesh for sanity check
%[analysis{2}.vertices,analysis{2}.faces] = create_hedgehog_mesh(analysis{2}.N, 0, 0.7);   %create "noisy sphere" mesh 

analysis{2}.vertices(:,1) = analysis{2}.vertices(:,1) + 2.5;      %shift second mesh along x-axis

hold off;
for i = 1:length(analysis);                %generate two meshes and corresponding algorithms
    analysis{i}.mesh = geodesic_new_mesh(analysis{i}.vertices, analysis{i}.faces);         %initilize new mesh
    analysis{i}.algorithm{1} = geodesic_new_algorithm(analysis{i}.mesh, 'exact');   %initialize exact algorithm
    analysis{i}.algorithm{2} = geodesic_new_algorithm(analysis{i}.mesh, 'subdivision', 4);   %initialize subdivision algorithm with 3 points per edge
    analysis{i}.algorithm{3} = geodesic_new_algorithm(analysis{i}.mesh, 'dijkstra');   %initialize dijkstra algorithm 
    
    trisurf(analysis{i}.faces,analysis{i}.vertices(:,1),analysis{i}.vertices(:,2),analysis{i}.vertices(:,3),...
        'FaceColor', 'w', 'EdgeColor', 'k', 'FaceAlpha', 0.99);       %plot the mesh
    hold on;
end;

colors = {'r','b', 'g'};           
for i = 1:length(analysis);
    vertex_id = 1;                      %put a single source at vertex #1
    source_points = {geodesic_create_surface_point('vertex',vertex_id, analysis{i}.vertices(vertex_id,:))};
    
    vertex_id = analysis{i}.N;           %last vertex of the mesh is destination
    destination = geodesic_create_surface_point('vertex',vertex_id, analysis{i}.vertices(vertex_id,:));
    for j=1:length(analysis{i}.algorithm)
        current_algorithm = analysis{i}.algorithm{j};
        geodesic_propagate(current_algorithm, source_points);   %propagation stage of the algorithm 
        path = geodesic_trace_back(current_algorithm, destination); %find a shortest path from source to destination
        
        [x,y,z] = extract_coordinates_from_path(path);
        plot3(x,y,z,colors{j},'LineWidth',2);    %plot a sinlge path for this algorithm

        path_length = sum(sqrt(diff(x).^2 + diff(y).^2 + diff(z).^2));            %length of the path
        [source_id, d] = geodesic_distance_and_source(current_algorithm, destination);    %you can confirm that the path length is equal to the distance estimated by geodesic_best_source
        s = sprintf('mesh %d, %s algorithm, estimated/actual length of the path is %f/%f', ...
            i, current_algorithm.type, d, path_length);
        disp(s);
    end
end
legend('mesh 1','mesh 2', 'exact', 'subdivision', 'dijkstra');

xlim([-1 3.5]);
ylim([-1 1]);
zlim([-1 1]);
daspect([1 1 1]);

geodesic_delete;