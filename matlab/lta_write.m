function err = lta_write(fname,lta)
% err = lta_write(fname,ltastruct)
%
% This creates a legal LTA file given an LTA "struct". The elements of
% the LTA struct are below. Note that it is different than the LTA
% struct in the C code (ie, in include/transform.h). nxforms is
% always 1. mean and sigma are always 0.
% 
% lta.type  - 0=vox2vox
% lta.xform - the actual 4x4 matrix
% lta.srcfile - the source volume file name
% lta.srcmri  - the mri struct as read by MRIread(). The .vol
%   element is not necessary.
% The items below are the same but for the destination volume
%   lta.dstfile 
%   lta.dstmri
% lta.subject - subjectname (or just use '' if not known)
%

err = 1;
fp = fopen(fname,'w');
if(fp == -1)
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

fprintf(fp,'# transform file %s\n',fname);
fprintf(fp,'# created by lta_write.m\n');
fprintf(fp,'type     = %d\n',lta.type);
fprintf(fp,'nxforms  = %d\n',1);
fprintf(fp,'mean     = 0 0 0\n');
fprintf(fp,'sigma    = 0 \n');
fprintf(fp,'1 4 4\n');
for r = 1:4
  for c = 1:4
    fprintf(fp,'%18.15f ',lta.xform(r,c));
  end
  fprintf(fp,'\n');
end
for n = 1:2
  if(n==1) 
    fprintf(fp,'src volume info\n'); 
    mri = lta.srcmri;
    mrifname = lta.srcfile;
  end
  if(n==2) 
    fprintf(fp,'dst volume info\n'); 
    mri = lta.dstmri;
    mrifname = lta.dstfile;
  end
  fprintf(fp,'valid = 1\n');
  fprintf(fp,'filename = %s\n',mrifname);
  fprintf(fp,'volume = %d %d %d\n',mri.volsize(1),mri.volsize(2),mri.volsize(3));
  fprintf(fp,'voxelsize = %18.15f %18.15f %18.15f\n',mri.volres(1),mri.volres(2),mri.volres(3));
  fprintf(fp,'xras = %18.15f %18.15f %18.15f\n',mri.x_r,mri.x_a,mri.x_s);
  fprintf(fp,'yras = %18.15f %18.15f %18.15f\n',mri.y_r,mri.y_a,mri.y_s);
  fprintf(fp,'zras = %18.15f %18.15f %18.15f\n',mri.z_r,mri.z_a,mri.z_s);
  fprintf(fp,'cras = %18.15f %18.15f %18.15f\n',mri.c_r,mri.c_a,mri.c_s);
end
if(isempty(lta.subject)) lta.subject = 'unknown'; end
fprintf(fp,'subject %s\n',lta.subject);

err = 0;
return

  
