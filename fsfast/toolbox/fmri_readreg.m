function [regmat, subject, inres, betres, intensity] = fmri_readreg(regname)
%
% [regmat subject inres betres intensity] = fmri_readreg(regname)  
%
% subject - name of subject
% inres   - in-plane resolution (mm)
% betres  - between-plane resolution (mm)
% ireg    - intensity value for register program
% regmat  - 4x4 matrix which converts a mm location in corronal space
%             to a mm location in functional space.
%
%


%
% fmri_readreg.m
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

regmat = [];

if(nargin ~= 1) 
  msg = 'USAGE: [regmat subject inres betres intensity]  = fmri_readreg(regname)';
  qoe(msg); error(msg);
end

regname = deblank(regname);
fid = fopen(regname);
if(fid == -1)
  msg = sprintf('Could not open  %s',regname);
  qoe(msg); error(msg);
end

[subject count] = fscanf(fid,'%s',1);
inres     = fscanf(fid,'%f',1);
betres    = fscanf(fid,'%f',1);
intensity = fscanf(fid,'%f',1);
regmat    = fscanf(fid,'%f',inf);
regmat    = reshape(regmat, [4 4])'; %'
fclose(fid);


return;
