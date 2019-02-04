function result = isInGeodesicROI(mesh_total,wholepath,face)

wholepath=wholepath(2:end);
verticesinfaces=mesh_total.faces(face,:);
b = intersect(verticesinfaces,wholepath);

resultat = size(b,2) < 1;
result = resultat(1);
