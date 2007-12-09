function result = isVertexInRadius(vertex, origin, radius)

vt = vertex - origin;

result = ( ((vt(1)^2 + vt(2)^2) <= (radius^2)) & ...
            ((vt(2)^2 + vt(3)^2) <= (radius^2)) & ...
            ((vt(1)^2 + vt(3)^2) <= (radius^2)) );
        
     