function mri = MRIread(fstring,headeronly)
% mri = MRIread(fstring,headeronly)
%
% Reads in a volume based on the fstring. fstring can be:
% 1. A stem, in which case the format and full file name is determined
%     by finding a file on disk called fstring.ext, where ext can be
%     either mgh, mgz, or bhdr
%  2. MGH file. Eg, f.mgh or f.mgz
%  3. BVolume HDR file. Eg, f.bhdr 
%
% Creates a structure similar to the FreeSurfer MRI struct
% defined in mri.h. Times are in ms and angles are in radians.
% The vox2ras0 matrix is the matrix that converts a 0-based
% column, row, and slice to XYZ. vox2ras1 is the same with
% 1-based indices. The volume is rows, cols, slices frames,
% but the vox2ras expects col, row, slice.
%
% If headeronly=1, then the pixel data is not read in.
%
% $Id: MRIread.m,v 1.4 2004/11/13 16:44:28 greve Exp $

mri = [];

if(nargin < 1 | nargin > 2)
  fprintf('mri = MRIread(fstring,headeronly)\n');
  return;
end
if(exist('headeronly')~=1) headeronly = 0; end

[fspec fstem fmt] = MRIfspec(fstring);
if(isempty(fspec))
  fprintf('ERROR: cannot determine format of %s\n',fstring);
  return;
end

%-------------- MGH ------------------------%
switch(fmt)
  case {'mgh','mgz'}
  [mri.vol, M, mr_parms, volsz] = load_mgh(fspec,headeronly);
  if(isempty(M))
    fprintf('ERROR: loading %s as MGH\n',fspec);
    mri = [];
    return;
  end
  if(~headeronly)
    mri.vol = permute(mri.vol,[2 1 3 4]);
    volsz = size(mri.vol);
  else
    mri.vol = [];
    volsz(1:2) = [volsz(2) volsz(1)];
  end
  tr = mr_parms(1);
  flip_angle = mr_parms(2);
  te = mr_parms(3);
  ti = mr_parms(4);
%--------------- bshort/bfloat -----------------------%  
 case {'bhdr'}
  if(~headeronly)
    [mri.vol bmri] = fast_ldbslice(fstem);  
    if(isempty(mri.vol))
      fprintf('ERROR: loading %s as bvolume\n',fspec);
      mri = [];
      return;
    end
    volsz = size(mri.vol);
  else
    mri.vol = [];
    bmri = fast_ldbhdr(fstem);
    if(isempty(bmri))
      fprintf('ERROR: loading %s as bvolume\n',fspec);
      mri = [];
      return;
    end
    [nslices nrows ncols ntp] = fmri_bvoldim(fstem);
    volsz = [nrows ncols nslices ntp];
  end
  M = bmri.T;
  tr = bmri.tr;
  flip_angle = bmri.flip_angle;
  te = bmri.te;
  ti = bmri.ti;
 otherwise
  fprintf('ERROR: format %s not supported\n',fmt);
  mri = [];
  return;
end
%--------------------------------------%

mri.fspec = fspec;
mri.pwd = pwd;
mri.flip_angle = flip_angle;
mri.tr  = tr;
mri.te  = te;
mri.ti  = ti;

mri.vox2ras0 = M;
mri.vox2ras1 = vox2ras_0to1(M);

% Dimensions not redundant when using header only
volsz(length(volsz)+1:4) = 1; % Make sure all dims are represented
mri.volsize = volsz(1:3); % only spatial components
mri.height  = volsz(1);   % Note that height (rows) is the first dimension
mri.width   = volsz(2);   % Note that width (cols) is the second dimension
mri.depth   = volsz(3);
mri.nframes = volsz(4);

%--------------------------------------------------------------------%
% Everything below is redundant in that they can be derivied from
% stuff above, but they are in the MRI struct defined in mri.h, so I
% thought I would add them here for completeness.

mri.xsize = sqrt(sum(M(:,1).^2));
mri.ysize = sqrt(sum(M(:,2).^2));
mri.zsize = sqrt(sum(M(:,3).^2));
mri.volres = [mri.xsize mri.ysize mri.zsize];

mri.x_r = M(1,1)/mri.xsize;
mri.x_a = M(2,1)/mri.xsize;
mri.x_s = M(3,1)/mri.xsize;

mri.y_r = M(1,2)/mri.ysize;
mri.y_a = M(2,2)/mri.ysize;
mri.y_s = M(3,2)/mri.ysize;

mri.z_r = M(1,3)/mri.zsize;
mri.z_a = M(2,3)/mri.zsize;
mri.z_s = M(3,3)/mri.zsize;

% Matrix of direction cosines
mri.Mdc = [M(1:3,1)/mri.xsize M(1:3,2)/mri.ysize M(1:3,3)/mri.zsize];

ic = [(mri.width)/2 (mri.height)/2 (mri.depth)/2 1]';
c = M*ic;
mri.c_r = c(1);
mri.c_a = c(2);
mri.c_s = c(3);
%--------------------------------------------------%


return;











