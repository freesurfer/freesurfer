function [fspec, fstem, fmt] = MRIfspec(fstring)
% [fspec fstem fmt] = MRIfspec(fstring)
%
% Determine the file specification (fspec), stem (fstem), and format
% (fmt) from the given file string.
%
% A specification is the name of a file as it could exist on disk. The
% specification has the form fstem.fmt, where fmt can be mgh, mgz,
% bhdr, or img. fstring can be either an fspec or fstem. 
%
% If fstring is an fspec, then the format is determined from the
% extension.
%
% If fstring is an fstem, then the format is determined by finding a
% file on disk named fstem.fmt. The formats are searched in the
% following order: mgh, mgz, bhdr, or img. The search order is only
% important when there are multiple files that would have met the
% criteria, then only the first one is chosen. If no such file is
% found, then empty strings are returned.
%
% $Id: MRIfspec.m,v 1.1 2004/11/13 16:48:43 greve Exp $

fspec = [];
fstem = [];
fmt   = [];

if(nargin ~= 1)
  fprintf('[fspec fstem fmt] = MRIfspec(fstring)\n');
  return;
end

% First, examin fstring to see if it has an extension fspec must have
% at least 4 characters (5 for bhdr). Order is not important here.
if(length(fstring) > 4) 
  ext = fstring(end-2:end);
  switch(ext)
   case 'mgh',
    fspec = fstring;
    fmt = 'mgh';
    fstem = fstring(1:end-4);
    return;
   case 'mgz',
    fspec = fstring;
    fmt = 'mgz';
    fstem = fstring(1:end-4);
    return;
   case 'img',
    fspec = fstring;
    fmt = 'img';
    fstem = fstring(1:end-4);
    return;
  end
  if(length(fstring) > 5) 
    ext = fstring(end-3:end);
    switch(ext)
     case 'bhdr',
      fspec = fstring;
      fmt = 'bhdr';
      fstem = fstring(1:end-5);
      return;
    end
  end
end

% If it gets here, then it cannot determine the format from an
% extension, so fstring could be a stem, so see if a file with
% stem.fmt exists. Order is imporant here.
fstem = fstring;

fmt = 'mgh';
fspec = sprintf('%s.%s',fstring,fmt);
if(fast_fileexists(fspec)) return; end

fmt = 'mgz';
fspec = sprintf('%s.%s',fstring,fmt);
if(fast_fileexists(fspec)) return; end

fmt = 'bhdr';
fspec = sprintf('%s.%s',fstring,fmt);
if(fast_fileexists(fspec)) return; end

fmt = 'img';
fspec = sprintf('%s.%s',fstring,fmt);
if(fast_fileexists(fspec)) return; end

% If it gets here, then could not determine format
fstem = [];
fmt = [];
fspec = [];

return;