function s = ldstruct(fid)
% s = ldstruct(fid)
%
% Loads a structure saved with svstruct().  fid is a file
% identifier opened for reading.
%
%  $Id: ldstruct.m,v 1.1 2003/03/04 20:47:41 greve Exp $

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
     nstrings = vsize(1);
     v = [];
     for m = 1:nstrings,
       sline = fgetl(fid);
       v = strvcat(v,sline);  
     end
  end

  v = reshape(v, vsize);
  s = setfield(s,fld,v);

end

return
