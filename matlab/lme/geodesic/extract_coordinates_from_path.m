function [x,y,z] = extract_coordinates_from_path(path)

%this elegant way is incompatible with the older versions of matlab
% x = cellfun(@(p) p.x, path);    %the simplest way to extract coordinates from the path
% y = cellfun(@(p) p.y, path);    %if it looks complicated, you can use "for" similar to example1.m
% z = cellfun(@(p) p.z, path);

x = zeros(length(path),1);
y = x; 
z = y;

for i=1:length(path)
    x(i) = path{i}.x;
    y(i) = path{i}.y;
    z(i) = path{i}.z;
end;

