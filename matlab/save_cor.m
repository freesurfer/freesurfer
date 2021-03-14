function err = save_cor(corvol,stemdir,dir)
%
% err = save_cor(corvol,<stemdir>,<dir>)
%
% Saves in stemdir/dir. If stemdir and dir are not 
% specified, saves in the current directory.
%
% See also load_cor.
%


%
% save_cor.m
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


err = 1;

if(nargin < 1 | nargin > 3)
  msg = 'USAGE: err = save_cor(corvol,<stemdir>,<dir>)';
  qoe(msg);error(msg);
end

if(exist('stemdir') ~= 1 & exist('dir') ~= 1)
  cordir = '.';
else
  cordir = sprintf('%s/%s',stemdir,dir);
  status = mkdir(stemdir,dir);
end

fname = sprintf('%s/COR-.info',cordir);
fid = fopen(fname,'w');
if(fid == -1)
  fprintf('ERROR: opening %s for writing\n',fname);
  return;
end

fprintf(fid,'imnr0 %d\n',1);
fprintf(fid,'imnr1 %d\n',size(corvol,3));
fprintf(fid,'ptype %d\n',2);
fprintf(fid,'x %d\n',size(corvol,2));
fprintf(fid,'y %d\n',size(corvol,1));
fprintf(fid,'fov %f\n',0.256);
fprintf(fid,'thick %f\n',0.001);
fprintf(fid,'psiz %f\n',0.001);
fprintf(fid,'locatn %f\n',0);
fprintf(fid,'strtx %f\n',-0.128);
fprintf(fid,'endx %f\n',0.128);
fprintf(fid,'strty %f\n',-0.128);
fprintf(fid,'endy %f\n',0.128);
fprintf(fid,'strtz %f\n',-0.128);
fprintf(fid,'endz %f\n',0.128);
fprintf(fid,'tr %f\n',0.0);
fprintf(fid,'te %f\n',0.0);
fprintf(fid,'ti %f\n',0.0);
fclose(fid);

[nc,nr,ns] = size(corvol);

for s=1:ns
  corslice = squeeze(corvol(:,:,s))'; %' Convert to column major 
  corslice = uint8(max(0,min(255,corslice)));
  corslicefile = sprintf('%s/COR-%03d',cordir,s);
  fid=fopen(corslicefile,'wb');
  if(fid == -1)
    msg = sprintf('Could not open %s for writing.',corslicefile); 
    qoe(msg); error(msg);
  end
  precision = 'uint8';
  Nv = prod(size(corslice));
  count = fwrite(fid,corslice(:),precision);
  fclose(fid); 
  if(count ~= Nv)
    fprintf(2,'ERROR: wrote %d/%d elements to %s\n',count,Nv,corslicefile);
    err = 1;
    return;
  else 
    err = 0;
  end
end

err = 0;

return;
