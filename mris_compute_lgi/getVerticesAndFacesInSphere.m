function [verticesList, facesList] = getVerticesAndFacesInSphere(mesh_outer, iV, radius)

% find all the vertices included in a sphere
verticesList = [];

for vertex = 1:length(mesh_outer.vertices)
    if isVertexInRadius(mesh_outer.vertices(vertex,:), mesh_outer.vertices(iV,:), radius)
        verticesList = [verticesList,vertex];
    end
end

% find faces to which those vertices belong
facesList = [];
for vert = verticesList
    facesList = [facesList, mesh_outer.facesOfVertex(vert).faceList];
end

facesList = unique(facesList);
