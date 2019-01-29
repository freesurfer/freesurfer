function totalArea = getFacesArea(mesh_in, ListFaces)

% to measure the area of a ROI given by a list of faces

totalArea = 0;

for iF = 1:size(ListFaces,1)
    totalArea = totalArea + getFaceArea(mesh_in,ListFaces,iF);
end

