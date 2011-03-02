function mriimg = MRIextractImage(mrivol,index,dim)
% mriimg = MRIextractImage(mrivol,index,dim)
%
% Extracts an image from the volume within the dim at the
% given index and returns an mristruct with an accurate
% geometry. Returns all frames in mriimg as the 4th dim.
%
%  dim = 1 = row,   mriimg.vol = mrivol.vol(index,:,:);
%  dim = 2 = col,   mriimg.vol = mrivol.vol(:,index,:);
%  dim = 3 = slice, mriimg.vol = mrivol.vol(:,:,index);
%
% Eg, to extract an axial from orig.mgz:
% mrivol = MRIread('orig.mgz');
% mriimg = MRIextractImage(mrivol,128,1);
% MRIwrite(mriimg,'orig.slice.mgz');
% tkregister2 --mov orig.slice.mgz --targ orig.mgz \
%     --regheader --reg tmp.reg
% 
% Not 100% sure it works (very close though).
%


%
% MRIextractImage.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
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


mriimg = [];
if(nargin < 3 | nargin > 3)
  fprintf('mriimg = MRIextractImage(mrivol,index,dim)\n');
  return;
end

if(dim < 1 | dim > 3)
  fprintf('ERROR: dim = %d out of range\n',dim);
  return;
end
if(index < 1 | index > mrivol.volsize(dim))
  fprintf('ERROR: index = %d out of range.\n',index);
  fprintf('       range is 1 to \n',mrivol.volsize(dim));
  return;
end

mriimg = mrivol;
mriimg.vol = [];

crs0 = zeros(4,1);
crs0(4) = 1;
if(dim == 1)
  % row-image: made of cols and slices 
  % cols of img are vol slices
  % rows of img are vol cols
  mriimg.vol = squeeze(mrivol.vol(index,:,:,:));
  crs0(2) = index-1;
  mriimg.Mdc = mrivol.Mdc(:,[3 1 2]);
  mriimg.xsize = mrivol.zsize;
  mriimg.ysize = mrivol.xsize;
  mriimg.zsize = mrivol.ysize;
elseif(dim == 2)
  % col-image: made of rows and slices 
  % cols of img are vol slices
  % rows of img are vol rows
  mriimg.vol = squeeze(mrivol.vol(:,index,:,:));
  crs0(1) = index-1;
  mriimg.Mdc = mrivol.Mdc(:,[3 2 1]);
  mriimg.xsize = mrivol.zsize;
  mriimg.ysize = mrivol.ysize;
  mriimg.zsize = mrivol.xsize;
else
  % slice-image: made of rows and cols
  % cols of img are vol cols
  % rows of img are vol rows
  mriimg.vol = squeeze(mrivol.vol(:,:,index,:));
  % Can keep Mdc and {x,y,z}size copied from vol
  crs0(3) = index-1;
end

% Make sure frames are in 4th dim and 3rd dimsize=1
mriimg.vol = permute(mriimg.vol,[1 2 4 3]);

P0 = mrivol.vox2ras0(1:3,:) * crs0;

mriimg.volsize = size(mriimg.vol);
mriimg.volsize(3) = 1;
mriimg.volsize(4) = mrivol.nframes;
mriimg.height  = mriimg.volsize(1);
mriimg.width   = mriimg.volsize(2);
mriimg.depth   = mriimg.volsize(3);
mriimg.nframes = mriimg.volsize(4);
mriimg.volres  = [mriimg.xsize mriimg.ysize mriimg.zsize];
mriimg.vox2ras0 = [mriimg.Mdc*diag(mriimg.volres) P0; 0 0 0 1];
mriimg.vox2ras  = mriimg.vox2ras0;
mriimg.tkrvox2ras = vox2ras_tkreg(mriimg.volsize,mriimg.volres);
mriimg.vox2ras1 = vox2ras_0to1(mriimg.vox2ras0); 

return;

