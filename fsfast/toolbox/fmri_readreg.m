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
% $Id: fmri_readreg.m,v 1.1 2003/03/04 20:47:40 greve Exp $

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
