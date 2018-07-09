function [mosimg, mos_n_rows, mos_n_cols] = mosaic_tb(img, mosaic_cols)
%MOSAIC_TB Turns 3-D data into 2-D mosaic format with given number of columns.
%   [MOSIMG, MOS_N_ROWS, MOS_N_COLS] = MOSAIC_TB(IMG, MOSAIC_COLS),
%   where IMG is the 3-D image stack and MOSAIC_COLS are the number of columns 
%   of the mosaic. MOSIMG is the 2-D mosaic-ed image with MOS_N_ROWS and
%   MOS_N_COLS.
%
%   If the given data is found inconsistent an empty array MOSIMG and
%   a value of -1 is returned for MOS_N_ROWS, and MOS_N_COLS.
%
%   Copyright (C) 2003  Thomas Benner  <thomas.benner@nmr.mgh.harvard.edu>

mosimg = [];
mos_n_rows = -1;
mos_n_cols = -1;

% we only accept 3-D arrays
if (ndims(img) ~= 3)
	mosimg = img;
	mos_n_rows = size(img, 1);
	mos_n_cols = size(img, 2);
	return;
end

% do nothing if mosaic_cols < 2
if (mosaic_cols<=1)
	mosimg = img;
	mos_n_rows = size(img, 1);
	mos_n_cols = size(img, 2);
	return;
end

% how many rows are needed for mosaic image?
mosaic_rows = ceil(size(img, 3) / mosaic_cols);

mos_n_rows = size(img, 1) * mosaic_rows;
mos_n_cols = size(img, 2) * mosaic_cols;

mosimg = zeros(mos_n_rows, mos_n_cols);

% reformat array (mosaic)
for slc = 0:size(img, 3)-1
  row_start = fix(slc/mosaic_cols)  * size(img, 1) + 1;
  col_start = mod(slc, mosaic_cols) * size(img, 2) + 1;
  mosimg(row_start:row_start+size(img, 1)-1, col_start:col_start+size(img, 2)-1) ...
    = img(:,:,slc+1);
end

return;
