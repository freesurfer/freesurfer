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
%    $Date: 2011/03/02 00:04:06 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
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
