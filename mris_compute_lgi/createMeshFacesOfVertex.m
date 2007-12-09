function mesh_out = createMeshFacesOfVertex(vertices, faces)

%% Original Author: Lucas Tamarit, lucas.tamarit@cisa.unige.ch

mesh_out.vertices = vertices;
mesh_out.faces = faces;

mesh_out.facesOfVertex = repmat(struct('faceList',[]),size(vertices,1),1);

for iF = 1:size(faces,1)
    if (mod(iF,5000) == 0)
        disp(['face ',num2str(iF),' / ',num2str(size(faces,1))])
    end

    for iV = 1:3
        mesh_out.facesOfVertex(faces(iF,iV)).faceList = ...
    union( mesh_out.facesOfVertex(faces(iF,iV)).faceList, iF);
    end
end
