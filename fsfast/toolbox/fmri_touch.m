function r = fmri_touch(fname);
% r = fmri_touch(fname);
% 
% Simple function to create a file called fname.  This is supposed
% to be something like the unix touch, but it has
% no effect if fname already exists.
%
% $Id: fmri_touch.m,v 1.1 2003/03/04 20:47:40 greve Exp $

if(nargin ~= 1) 
  msg = 'USAGE: r = fmri_touch(fname);';
  qoe(msg); error(msg);
end

fname = deblank(fname);

fid = fopen(fname,'a');
if(fid == -1) 
  msg = sprintf('Could not open %s for appending',fname);
  qoe(msg); error(msg);
end

fclose(fid);

r = 0;
return;
