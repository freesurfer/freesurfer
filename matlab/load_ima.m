function [imvol,M,h] = load_ima(fnamestem)
% [imvol,M,h] = load_ima(fnamestem)


%
% load_ima.m
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



filepaths = ls('-1v',sprintf('%s*.ima',fnamestem));

filelist = [];
rest = filepaths;
while ~isempty(rest)
  [tok,rest] = strtok(rest);
  if ~isempty(tok)
    filelist = strvcat(filelist,tok);
  end
end

nslices = size(filelist,1);

for i = [1:nslices]
  fname = deblank(filelist(i,:));
  h = read_siemens_header(fname);
  nrows = h.h_G28_Pre_Rows;
  ncols = h.h_G28_Pre_Columns;
  height = h.h_G21_Rel1_CM_FoV_Height;
  width = h.h_G21_Rel1_CM_FoV_Width;
  slice_thickness = h.h_G18_Acq_SliceThickness;
  pos_vec = [h.h_G21_Rel1_CM_ImagePosition_Sag, h.h_G21_Rel1_CM_ImagePosition_Cor, h.h_G21_Rel1_CM_ImagePosition_Tra];
  row_vec = [h.h_G21_Rel1_CM_ImageRow_Sag, h.h_G21_Rel1_CM_ImageRow_Cor, h.h_G21_Rel1_CM_ImageRow_Tra];
  col_vec = [h.h_G21_Rel1_CM_ImageColumn_Sag, h.h_G21_Rel1_CM_ImageColumn_Cor, h.h_G21_Rel1_CM_ImageColumn_Tra];
  norm_vec = [h.h_G21_Rel1_CM_ImageNormal_Sag, h.h_G21_Rel1_CM_ImageNormal_Cor, h.h_G21_Rel1_CM_ImageNormal_Tra];
  ManufacturerModel = h.h_G51_Txt_ManufacturerModel;
  if i==1
    imvol = zeros(nrows,ncols,nslices);
    pos_vec0 = pos_vec;
  end
  imvol(:,:,i) = read_siemens_image(fname)';
  i = i+1;
end

%slice_spacing = (1+h.h_G21_Rel2_Mr_CurrentSliceDista)*slice_thickness

pos_vec1 = pos_vec;
if (nslices>1)
  slice_spacing = (pos_vec1-pos_vec0)*norm_vec'/(nslices-1);
else
  slice_spacing = slice_thickness;
end

M = [height/nrows*col_vec(1), width/ncols*row_vec(1), slice_spacing*norm_vec(1), pos_vec0(1) + ...
       height/nrows*row_vec(1)*(-0.5-0.5*nrows) + width/ncols*col_vec(1)*(-0.5-0.5*ncols) - slice_spacing*norm_vec(1); 
     height/nrows*col_vec(2), width/ncols*row_vec(2), slice_spacing*norm_vec(2), pos_vec0(2) + ...
       height/nrows*row_vec(2)*(-0.5-0.5*nrows) + width/ncols*col_vec(2)*(-0.5-0.5*ncols) - slice_spacing*norm_vec(2); 
     height/nrows*col_vec(3), width/ncols*row_vec(3), slice_spacing*norm_vec(3), pos_vec0(3) + ...
       height/nrows*row_vec(3)*(-0.5-0.5*nrows) + width/ncols*col_vec(3)*(-0.5-0.5*ncols) - slice_spacing*norm_vec(3);
     0, 0, 0, 1];

