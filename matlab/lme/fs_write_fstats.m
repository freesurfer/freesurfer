function fs_write_fstats(fstats,mri,fname,data_type)
% fs_write_fstats(fstats,mri,fname,data_type)
% 
% Writes F-statistic or p-value or Freesurfer significance maps to Freesurfer's 
% .mgh or .mgz data file. This can be useful for visualization and post-processing 
% in Freesurfer. 
%
% Input
% fstats: Structure obtained from lme_mass_F.
% mri: Mri structure (read with fs_read_Y).
% fname: Output file name.
% data_type: Determines what is going to be written. This input can be one 
% of three strings: 'fval' (signed F-statistic map), 'pval' (signed p-value map) 
% or 'sig' (Freesurfer significance map -log10(pval).*sgn).
%
% Original Author: Jorge Luis Bernal Rusiel 
%
if nargin < 4
    error('Too few inputs');
end;
mri.volsz(4) = 1;
if strcmpi(data_type,'fval')
    fs_write_Y(fstats.F.*fstats.sgn,mri,fname);
elseif strcmpi(data_type,'pval')
    pval = fstats.pval.*fstats.sgn;
    pval(pval==1) = 1;
    fs_write_Y(pval,mri,fname);
elseif strcmpi(data_type,'sig')
    fs_write_Y(-log10(fstats.pval).*fstats.sgn,mri,fname);
else
    error('Valid strings for data_type are ''fval'' or ''pval'' or ''sig''');
end;


