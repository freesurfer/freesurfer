function [vol,M] = load_mgh(fname)
%
% Loads the indicated corronal slices from 
%   $SUBJECTS_DIR/subject/mri/volumeid
%
% [vol, Mvox_to_ras] = load_mgh(fname)                
%
% $Id: load_mgh.m,v 1.1 2002/10/30 00:22:27 fischl Exp $

if(nargin < 1 | nargin > 1)
  msg = 'USAGE: vol = load_mgh(fname)';
  qoe(msg);error(msg);
end

fid = fopen(fname, 'rb', 'b') ;
v = fread(fid, 1, 'int') ; 
width = fread(fid, 1, 'int') ; 
height = fread(fid, 1, 'int') ; 
depth = fread(fid, 1, 'int') ; 
nframes = fread(fid, 1, 'int') ; 
type = fread(fid, 1, 'int') ; 
dof = fread(fid, 1, 'int') ; 

UNUSED_SPACE_SIZE= 256;
USED_SPACE_SIZE = (3*4+4*3*4);  % space for ras transform

unused_space_size = UNUSED_SPACE_SIZE-2 ;
ras_good_flag = fread(fid, 1, 'short') ; 
if (ras_good_flag)
	 unused_space_size = unused_space_size - USED_SPACE_SIZE ;
	 xsize = fread(fid, 1, 'float32') ; 
	 ysize = fread(fid, 1, 'float32') ; 
	 zsize = fread(fid, 1, 'float32') ; 

	 x_r = fread(fid, 1, 'float32') ; 
	 x_a = fread(fid, 1, 'float32') ; 
	 x_s = fread(fid, 1, 'float32') ; 

	 y_r = fread(fid, 1, 'float32') ; 
	 y_a = fread(fid, 1, 'float32') ; 
	 y_s = fread(fid, 1, 'float32') ; 
		
	 z_r = fread(fid, 1, 'float32') ; 
	 z_a = fread(fid, 1, 'float32') ; 
	 z_s = fread(fid, 1, 'float32') ; 
		
	 c_r = fread(fid, 1, 'float32') ; 
	 c_a = fread(fid, 1, 'float32') ; 
	 c_s = fread(fid, 1, 'float32') ; 
	 M(1,1) = xsize * x_r ;
	 M(1,2) = ysize * y_r ;
	 M(1,3) = zsize * z_r ;

	 M(2,1) = xsize * x_a ;
	 M(2,2) = ysize * y_a ;
	 M(2,3) = zsize * z_a ;

	 M(3,1) = xsize * x_s ;
	 M(3,2) = ysize * y_s ;
	 M(3,3) = zsize * z_s ;

	 ci = (width-1)/2 ; cj = (height-1)/2 ; ck = (depth-1)/2 ;
	 M(1,4) = c_r - (M(1,1)*ci + M(1,2)*cj + M(1,3)*ck) ;
   M(2,4) = c_a - (M(2,1)*ci + M(2,2)*cj + M(2,3)*ck);
   M(3,4) = c_s - (M(3,1)*ci + M(3,2)*cj + M(3,3)*ck);
	 M(4,4) = 1 ;
end

MRI_UCHAR =  0 ;
MRI_INT =    1 ;
MRI_LONG =   2 ;
MRI_FLOAT =  3 ;
MRI_SHORT =  4 ;
MRI_BITMAP = 5 ;

switch type
			 case MRI_FLOAT,
			   bpv = 4 ;
			 case MRI_UCHAR,
			   bpv = 1 ;
			 case MRI_SHORT,
			   bpv = 2 ;
			 case MRI_INT,
			   bpv = 4 ;
end

						
nelts = width * height ;  % bytes per slice 

fseek(fid, unused_space_size, 'cof') ;
vol = zeros(width, height, depth) ;

for z=1:depth
		switch type
			 case MRI_FLOAT,
					slice = fread(fid, nelts, 'float32') ; 
			 case MRI_UCHAR,
					slice = fread(fid, nelts, 'uchar') ; 
			 case MRI_SHORT,
					slice = fread(fid, nelts, 'short') ; 
			 case MRI_INT,
					slice = fread(fid, nelts, 'int') ; 
		end
		vol(:,:,z) = reshape(slice, [width height])' ;
end

fclose(fid) ;


