function s = ldstruct(fid)
% s = ldstruct(fid)
%
% Loads a structure saved with svstruct().  fid is a file
% identifier opened for reading. This assumes a format like:
%  fieldname type dim size
%    val1
%    <val2> ...
%  type is char or nonchar
%


%
% ldstruct.m
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

s = [];

if(nargin ~= 1) 
  msg = 'USAGE: s = ldstruct(fid)';
  qoe(msg); error(msg);
end

nfields = fscanf(fid,'%*s %d',1);

for n = 1:nfields,
  fld   = fscanf(fid,'%s',1);
  vtype = fscanf(fid,'%d',1);
  vdim  = fscanf(fid,'%d',1);
  vsize = fscanf(fid,'%d ',vdim)'; %'
  nv = prod(vsize);

  if(vtype == 0)
     %% Not a character string %%
     v = fscanf(fid,'%f',nv);
  else
     %% Character string %%
     % Swallow new line. I added this on 8/11/06. It was not broken
     % before but matlab 7.2 seems not to work without it. So, I
     % don't know if this will break previous versions.
     % NJS: it does break previous versions, so this check for version
     % number has been added before eating newline.
     % Note: v7.3 does not work correctly when newline is eaten!
     matlab_version = version;
     matlab_version = str2num(matlab_version(1:3));
     if (matlab_version == 7.2)
       fgetl(fid);
     end
     nstrings = vsize(1);
     v = [];
     for m = 1:nstrings,
       sline = fgetl(fid);
       v = strvcat(v,sline);  
     end
  end
  
  if(prod(vsize) ~= prod(size(v)) )
    fprintf('ERROR: size mismatch in field %s, ',fld);
    fprintf('%d vs %d\n',prod(vsize),prod(size(v)));
    s = [];
    return;
  end
  
  v = reshape(v, vsize);
  s = setfield(s,fld,v);
end

return
