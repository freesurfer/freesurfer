function [] = save_mgh(vol,fname, M, sizes)
%
% Loads the indicated mgh format volume
%   $SUBJECTS_DIR/subject/mri/volumeid
%
% save_mgh(vol,fname, M, sizes)  
%
% $Id: save_mgh.m,v 1.1 2002/10/30 00:22:28 fischl Exp $

if(nargin < 2 | nargin > 4)
  msg = 'USAGE: save_mgh(vol,fname)';
  qoe(msg);error(msg);
end

MRI_UCHAR =  0 ;
MRI_INT =    1 ;
MRI_LONG =   2 ;
MRI_FLOAT =  3 ;
MRI_SHORT =  4 ;
MRI_BITMAP = 5 ;
MRI_TENSOR = 6 ;
slices   = [1:256];

fid = fopen(fname, 'wb', 'b') ;
[width,height,depth,rows,cols] = size(vol) ;
fwrite(fid, 1, 'int') ;					 % magic #
fwrite(fid, width, 'int') ; 
fwrite(fid, height, 'int') ; 
fwrite(fid, depth, 'int') ; 
fwrite(fid, 1, 'int') ;					 % # of frames
if (ndims(vol) == 5)
	 is_tensor = 1 ;
	 fwrite(fid, MRI_TENSOR, 'int') ;	 % type = MRI_TENSOR
else
	 is_tensor = 0 ;
	 fwrite(fid, MRI_FLOAT, 'int') ;	 % type = MRI_FLOAT
end

fwrite(fid, 1, 'int') ;          % dof (not used)
dof = fread(fid, 1, 'int') ; 

UNUSED_SPACE_SIZE= 256;
USED_SPACE_SIZE = (3*4+4*3*4);  % space for ras transform

unused_space_size = UNUSED_SPACE_SIZE-2 ;
if (nargin > 2)
	 fwrite(fid, 1, 'short') ;       % ras_good_flag = 0
	 unused_space_size = unused_space_size - USED_SPACE_SIZE ;
	 fwrite(fid, sizes(1), 'float32') ; % xsize
	 fwrite(fid, sizes(2), 'float32') ; % ysize
	 fwrite(fid, sizes(3), 'float32') ; % zsize
	 
	 fwrite(fid, M(1,1), 'float32') ;  % x_r
	 fwrite(fid, M(2,1), 'float32') ;  % x_a
	 fwrite(fid, M(3,1), 'float32') ;  % x_s

	 fwrite(fid, M(1,2), 'float32') ;  % y_r
	 fwrite(fid, M(2,2), 'float32') ;  % y_a
	 fwrite(fid, M(3,2), 'float32') ;  % y_s

	 fwrite(fid, M(1,3), 'float32') ;  % z_r
	 fwrite(fid, M(2,3), 'float32') ;  % z_a
	 fwrite(fid, M(3,3), 'float32') ;  % z_s

	 fwrite(fid, M(1,4), 'float32') ;  % c_r
	 fwrite(fid, M(2,4), 'float32') ;  % c_a
	 fwrite(fid, M(3,4), 'float32') ;  % c_s
else
	 fwrite(fid, 0, 'short') ;       % ras_good_flag = 0

end

fwrite(fid, zeros(unused_space_size,1), 'char') ;
bpv = 4 ;   % bytes/voxel

						
nelts = width * height ;  % bytes per slice 

for row=1:rows
	for col=1:cols
		for z=1:depth
				slice = reshape(vol(:,:,z, row, col), [width height]) ;
				fwrite(fid, slice, 'float32') ; 
				vol(:,:,z) = reshape(slice, [width height]) ;
		end
	end
end

fclose(fid) ;


