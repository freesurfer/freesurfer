function err = write_wfile(fname, w, v)
% err = write_wfile(fname, w, <v>)
% 
% writes a vector into a binary 'w' file
%  fname - name of file to write to
%  w     - vector of values to be written
%  v     - 0-based vertex numbers 
%          (assumes 0 to N-1 if not present or empty).
%
% See also read_wfile.
%


%
% write_wfile.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:13 $
%    $Revision: 1.3 $
%
% Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%


err = 1;

if(nargin ~= 2 & nargin ~= 3)
  fprintf('USAGE: err = write_wfile(fname, w, <v>) \n');
  return;
end

vnum = length(w) ;

% Handle when v not given or is empty %
if(exist('v') ~= 1) v = []; end
if(isempty(v)) v = [0:vnum-1]; end

% open it as a big-endian file
fid = fopen(fname, 'wb', 'b') ;
if(fid == -1)
  fprintf('ERROR: could not open %s\n',fname);
  return;
end

fwrite(fid, 0, 'int16') ;
fwrite3(fid, vnum) ;
for i=1:vnum
  fwrite3(fid, v(i)) ;          % vertex number (0-based)
  fwrite(fid,  w(i), 'float') ; % vertex value
end

fclose(fid) ;

err = 0;

return
