

%
% read_label_old.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:10 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
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
