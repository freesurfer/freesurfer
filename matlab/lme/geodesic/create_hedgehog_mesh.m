%"smoothness" should be specified between 0(smooth, convex mesh) and 1 (a lot of sharp features)
%"waist" should be between 0 and 1; 0 means spherical mesh
function [p,tri] = create_hedgehog_mesh(N, smoothness, waist)

p = rand(N,3) - 0.5;
for i=1:N;
    p(i,:) = p(i,:)/norm(p(i,:));
end;

tri = convhulln(p);     

if nargin < 2
    smoothness = 0;
end

if nargin < 3
    waist = 0;
end

for i=1:N;
    scale = 1 + smoothness*(0.5-rand);
    p(i,:) = p(i,:)*scale;                      % add radial noise
    if max(abs(p(i,:))) > 1
        p(i,:) = p(i,:)/max(abs(p(i,:)));
    end
    
    scale = abs(p(i,1))*waist + (1 - waist);
    p(i,2:3) = p(i,2:3)*scale;                      % add "waist" on x-axis
end;
