

%
% read_label_old.m
%
% Original Author: Bruce Fischl
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%



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
