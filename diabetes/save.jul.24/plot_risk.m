g = 1:500;

a = .2370 ;
b = -36.21 ;
c = 6.0e-5;
d = 177 ;

r = a*g ;
ind = find(g<d);
r(ind) = r(ind) + b + (c * (d - g(ind)).^3);
r = a*g + b + c*((g<=d).*((d-g).^3));
plot(g,r);
