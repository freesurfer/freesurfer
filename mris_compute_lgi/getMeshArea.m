function area = getMeshArea (mesh)

area=0;
for i=1:size(mesh.faces)
    p1=mesh.faces(i,1);
    p2=mesh.faces(i,2);
    p3=mesh.faces(i,3);
    first_comp=[mesh.vertices(p1,2) mesh.vertices(p1,3) 1 ; mesh.vertices(p2,2) mesh.vertices(p2,3) 1; mesh.vertices(p3,2) mesh.vertices(p3,3) 1];
    scnd_comp=[mesh.vertices(p1,3) mesh.vertices(p1,1) 1 ; mesh.vertices(p2,3) mesh.vertices(p2,1) 1; mesh.vertices(p3,3) mesh.vertices(p3,1) 1];
    third_comp=[mesh.vertices(p1,1) mesh.vertices(p1,2) 1 ; mesh.vertices(p2,1) mesh.vertices(p2,2) 1; mesh.vertices(p3,1) mesh.vertices(p3,2) 1];
    a=det(first_comp);
    a=a^2;
    b=det(scnd_comp);
    b=b^2;
    c=det(third_comp);
    c=c^2;
    area_i=0.5*(sqrt(a+b+c));
    area=area+area_i;
end

