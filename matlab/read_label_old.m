%fname_label='lh-v1l.label';
fid=fopen(fname_label);
s = fgetl(fid);
s = fgetl(fid);
[num_vert_label] = sscanf(s,'%d');
vert_num = zeros(num_vert_label,1);
vert_data = zeros(5);
for vert = 1:1:num_vert_label,
  s = fgetl(fid);
  vert_data = sscanf(s,'%d %f %f %f %f');
  vert_num(vert) = vert_data(1:1);
end;
fclose(fid);
