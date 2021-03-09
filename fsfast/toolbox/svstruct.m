function svstruct(strct, fid)
% svstruct(strct, fid)
%
% Saves a structure as text to file.  fid is a file id opened for
% writing. The format is:
%
%    fieldname type dim size 
%    value(s)
%
%  type = 0 for numeric, 1 for char
%  dim is the number of dimensions in the data field (<3)
%  size is a dim-length vector of the number of elements 
%     in each dimension.
%
%  See also: ldstruct()
%


%
% svstruct.m
%
% Original Author: Doug Greve
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

if(nargin ~= 2) 
  msg = 'USAGE: svstruct(strct,fid)';
  qoe(msg); error(msg);
end

flds = fieldnames(strct);
fprintf(fid,'nfields %d\n',length(flds));
for n = 1:length(flds);
  fld = flds{n};
  v     = getfield(strct,fld);
  vsize = size(v);
  vdim  = length(vsize);
  if(vdim > 2) 
     msg = sprintf('Field %s has dimension %d, too large',fld,vdim);
     fclose(fid);
  end
  vtype = ischar(v);
  fprintf(fid,'%s %d %d ',fld,vtype,vdim);
  fprintf(fid,'%d ',vsize);
  fprintf(fid,'\n',vsize);

  if(vtype == 0)
     %% Not a character string %%
     fprintf(fid,'%g ',v);
     fprintf(fid,'\n');
  else
     %% Character string %%
     nstrings = vsize(1);
     for m = 1:nstrings,
       fprintf(fid,'%s\n',v(m,:));
     end

  end

end

return
