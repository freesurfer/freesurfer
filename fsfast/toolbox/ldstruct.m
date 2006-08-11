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
%  $Id: ldstruct.m,v 1.2 2006/08/11 18:29:54 greve Exp $

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
     fgetl(fid);
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
