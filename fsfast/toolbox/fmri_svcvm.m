function fmri_svcvm(cvmstruct,stem)
% fmri_svcvm(cvmstruct,stem)


%
% fmri_svcvm.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:02:33 $
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

if(nargin ~= 2)
     msg = 'USAGE: fmri_svcvm(cvmstruct,stem)';
     qoe(msg); error(msg);
end

cvmstruct.sz = size(cvmstruct.cvm,1);

fname = sprintf('%s.bfloat',stem);
fmri_svbfile(cvmstruct.cvm,fname);
fname = sprintf('%s.cvm',stem);
fid = fopen(fname,'w');
if(fid == -1) 
  msg = sprintf('Could not open %s for writing',fname);
  qoe(msg);error(msg);
end
fprintf(fid,'%d\n',cvmstruct.version);
fprintf(fid,'%d\n',cvmstruct.n);
fprintf(fid,'%f\n',cvmstruct.d);
fprintf(fid,'%d\n',cvmstruct.sz);
fprintf(fid,'%d\n',cvmstruct.norm);
fprintf(fid,'%d\n',cvmstruct.inv);
if(~isempty(cvmstruct.acoravg) & ~isempty(cvmstruct.acorstd))  
  fprintf(fid,'%f\n',cvmstruct.acoravg);
  fprintf(fid,'%f\n',cvmstruct.acorstd);
end

fclose(fid);

return;
