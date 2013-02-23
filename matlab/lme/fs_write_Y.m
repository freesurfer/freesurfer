function fs_write_Y(Y,mri,fname)
% fs_write_Y(Y,mri,fname)
% 
% Writes volume Y to a Freesurfer's .mgh or .mgz data file.  
%
% Input
% Y: Data matrix or vector.
% mri: Mri structure (read with fs_read_Y).
% fname: Output file name.
%
% $Revision: 1.1 $  $Date: 2013/02/23 21:05:16 $
% Original Author: Jorge Luis Bernal Rusiel 
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2013/02/23 21:05:16 $
%    $Revision: 1.1 $
%
if nargin < 3
    error('Too few inputs');
end;
save_mgh(reshape(Y',mri.volsz),fname,mri.M,mri.mr_parms);