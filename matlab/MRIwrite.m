function err = MRIwrite(mri,fspec)
% err = MRIwrite(mri,fspec)
%
% Writes an MRI volume based on the fspec. fspec can be:
%  1. MGH file. Eg, f.mgh or f.mgz
%  2. BVolume Stem. Eg, f for f_000.bshort or f_000.bfloat
%
% mri should be a structure like that read by MRIread.m
%
% $Id: MRIwrite.m,v 1.1 2004/11/09 19:16:56 greve Exp $

err = 1;

if(nargin ~= 2)
  fprintf('err = MRIwrite(mri,fspec)\n');
  return;
end

%-------------- MGH ------------------------%
if(MRIisMGH(fspec)) 
  M = mri.vox2ras0;
  mr_parms = [mri.tr mri.flip_angle mri.te mri.ti];
  err = save_mgh(permute(mri.vol,[2 1 3 4]), fspec, M, mr_parms);  
 return;
end

%-------------- bfloat ------------------------%
bmri.te = mri.te;
bmri.tr = mri.tr;
bmri.ti = mri.ti;
bmri.flip_angle = mri.flip_angle;
bmri.voldim = [size(mri.vol,1) size(mri.vol,2) size(mri.vol,3)];
bmri.nframes = size(mri.vol,4);
bmri.T = mri.vox2ras0;
bmri.volres = [mri.xsize mri.ysize mri.zsize];
bmri.cdc = [mri.x_r mri.x_a mri.x_s];
bmri.rdc = [mri.y_r mri.y_a mri.y_s];
bmri.sdc = [mri.z_r mri.z_a mri.z_s];
bmri.P0 = mri.vox2ras0(4,1:3);
bmri.c  = [mri.c_r mri.c_a mri.c_s];
err = fast_svbslice(mri.vol,fspec,[],'',bmri);

return;











