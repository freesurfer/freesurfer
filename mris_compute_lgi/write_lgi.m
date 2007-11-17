function ok = write_lgi (mesh_pial, lgi, outputfile)
% ok = write_lGI (mesh_pial, lgi, outputfile)
%
%
% modified from write_label.m by Marie Schaer
%
% Example: write_lgi('lh.pial',lgi,'lh.lgi.asc')
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/11/17 01:17:32 $
%    $Revision: 1.1 $
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

ok = 0;

lxyz=mesh_pial.vertices;

if(~isempty(lgi) & ~isempty(lxyz))
  npoints1 = length(lgi); 
  npoints2 = size(lxyz,1); 
  if(npoints1 ~= npoints2)
    fprintf('ERROR: lGI file does not have the same number of vertices than the pial surface.\n');
    return;
  end
  npoints = npoints1;
elseif(~isempty(lgi))
  npoints = length(lgi); 
  lxyz = zeros(npoints,3); 
elseif(~isempty(lxyz))
  npoints = length(lxyz); 
  lgi = zeros(npoints,1); 
end

% open as an ascii file
fid = fopen(outputfile, 'w') ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',labelfile);
  return;
end

num = [1:npoints]; 
num = num';
lgi = reshape(lgi,[npoints 1]);
lxyz   = reshape(lxyz,[npoints 3]);

l = [num lxyz lgi];
fprintf(fid,'%u %f %f %f %f \n',l');
fprintf(fid,'ENDPATH');

fclose(fid) ;

ok = 1;

return;
