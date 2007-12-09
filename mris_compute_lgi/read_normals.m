function norm = read_normals (mesh_in, fname)

vnum = length(mesh_in.vertices);
fnum = length(mesh_in.faces);
tnum = vnum + fnum;

fid=fopen(fname,'r');
normals = textscan (fid, '%f %f %f %f',tnum,'headerlines',2);
fclose(fid);

norm = [normals{1} normals{2} normals{3}]; 
 
norm = norm(1:vnum,:);