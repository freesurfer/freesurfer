function vo = getOrthogonalVector(v)

vn = v ./ norm(v);

vo = [0,0,0];
while (isequal(vo,[0,0,0]))    
    vo = rand(1,3)-0.5;
    vc = vo*vn';
    vo = vo-(vn*vc);
end

vo = vo ./ norm(vo);
