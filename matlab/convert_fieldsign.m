function err = convert_fieldsign(infname,outfname)
% err = convert_fieldsign(infname,outfname) 
%
% Loads a field sign (infname) as saved by the void write_fieldsign()
% in tksurfer and saves in the given output format. This allows the
% field sign map to be loaded as an overlay in tksurfer. This can
% also be used to convert a field sign mask.
% 
% Eg:
%  In matlab:
%    convert_fieldsign('fieldsign-lh','lh.fieldsign.mgh');
%    convert_fieldsign('fieldsignmask-lh','lh.fieldsignmask.mgh');
%  In the unix shell: 
%    mri_mask lh.fieldsign.mgh lh.fieldsignmask.mgh lh.fieldsign-masked.mgh
%    tksurfer subjectname lh inflated \
%      -overlay lh.fieldsign-masked.mgh -fthresh 0.5
%
% $Id: convert_fieldsign.m,v 1.1 2008/05/13 22:22:37 greve Exp $

err = 1;
fp = fopen(infname,'r','ieee-be');
if(fp == -1)
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

vnum = fread(fp,1,'int32');
%fprintf('vnum = %d\n',vnum);

fs = fread(fp,inf,'float32');
fclose(fp);

if(length(fs) ~= vnum)
  fprintf('ERROR: number of data points (%d) does not equal vnum (%d)\n',...
	  length(fs),vnum);
  return;
end

mri.nframes = vnum; 
mri.vol = fs(:)';
mri.vox2ras0 = eye(4);
mri.volsize = [1 vnum 1];
mri.volres = [1 1 1];
mri.tr = 0;
mri.flip_angle = 0;
mri.te = 0;
mri.ti = 0;
mri.xsize = mri.volres(1);
mri.ysize = mri.volres(2);
mri.zsize = mri.volres(3);
MRIwrite(mri,outfname);

err = 0;

return