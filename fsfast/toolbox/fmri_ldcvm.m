function cvmstruct = fmri_ldcvm(stem,init)
% cvmstruct = fmri_ldcvm(stem, <init>)
% If init equals 1 then an empty structure is
% returned if the relevant files do not exist.
% Otherwise an error is generated.


%
% fmri_ldcvm.m
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

if(nargin ~= 1 & nargin ~= 2)
   msg = 'USAGE: cvmstruct = fmri_ldcvm(stem,<init>)';
   qoe(msg); error(msg);
end

if(nargin == 1) init = 0; end

cvmstruct = fmri_cvmstruct;

fname = sprintf('%s.cvm',stem);
fid = fopen(fname,'r');
if(fid == -1 & init == 1) return; end
if(fid == -1)
  msg = sprintf('Could not open %s for reading',fname);
  qoe(msg);error(msg);
end

cvmstruct.version = fscanf(fid,'%d',1);
cvmstruct.n = fscanf(fid,'%d',1);
cvmstruct.d = fscanf(fid,'%f',1);

if(cvmstruct.version == 2)
  cvmstruct.sz   = fscanf(fid,'%d',1);
  cvmstruct.norm = fscanf(fid,'%d',1);
  cvmstruct.inv  = fscanf(fid,'%d',1);
  if(cvmstruct.norm)
    cvmstruct.acoravg = fscanf(fid,'%f',cvmstruct.sz);
    cvmstruct.acorstd = fscanf(fid,'%f',cvmstruct.sz);
  end
end

fclose(fid);

fname = sprintf('%s.bfloat',stem);
cvmstruct.cvm = fmri_ldbfile(fname);

return;
