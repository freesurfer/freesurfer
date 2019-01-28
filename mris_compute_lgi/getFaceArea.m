function farea = getFaceArea(mesh_in,ListFaces, face)

% to measure the area of a single face (given by the index of a face in a list of faces)

v1 = ListFaces(face,1);
v2 = ListFaces(face,2);
v3 = ListFaces(face,3);

first_comp = [mesh_in.vertices(v1,2), mesh_in.vertices(v1,3), 1; ...
              mesh_in.vertices(v2,2), mesh_in.vertices(v2,3), 1; ...
              mesh_in.vertices(v3,2), mesh_in.vertices(v3,3) 1];
        
scnd_comp = [mesh_in.vertices(v1,3), mesh_in.vertices(v1,1), 1; ...
             mesh_in.vertices(v2,3), mesh_in.vertices(v2,1), 1; ...
             mesh_in.vertices(v3,3), mesh_in.vertices(v3,1), 1];
       
third_comp = [mesh_in.vertices(v1,1), mesh_in.vertices(v1,2), 1; ...
              mesh_in.vertices(v2,1), mesh_in.vertices(v2,2), 1; ...
              mesh_in.vertices(v3,1), mesh_in.vertices(v3,2), 1];
        
a = det(first_comp)^2;
b = det(scnd_comp)^2;
c = det(third_comp)^2;

farea = 0.5*(sqrt(a+b+c));
