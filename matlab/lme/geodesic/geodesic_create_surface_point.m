function p = geodesic_create_surface_point(type,id,x,y,z)

p.type = type;
p.id = id;

if nargin == 3 
    p.x = x(1);
    p.y = x(2);
    p.z = x(3);
else
    p.x = x;
    p.y = y;
    p.z = z;
end
