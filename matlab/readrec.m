function A = readrec(infile,type,resx,resy,slice,offset)
% READREC reads a rectangular slice from a raw infile and
% assigns it to a matrix.
%
%       a = READREC('infile','type',resx,resy,slice,offset):
%       infile = raw data set name must be
%                in single quotes.
%       type = data type must also be in
%              single quotes and can only
%              be: 'short'(2byt),'long'(4byt),
%                  'float'(4byt),'double'(8byt),
%                  'char'(1byt),'uchar'(1byt),
%                  'uint8' (unsigned 8 bit integer data).
%       resx  = pixels per line
%       resy  = lines
%       slice = slice number desired (assuming the rest of the 
%                               frames are of the same size)
%       offset = header offset in increments of the type defined
%
% E.g.
%       a = readrec('mammo.raw','uint8',986,987,1,0);
						     %

%
% readrec.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:10 $
%    $Revision: 1.2 $
%
% Copyright (C) 2002-2007,
% The General Hospital Corporation (Boston, MA). 
% All rights reserved.
%
% Distribution, usage and copying of this software is covered under the
% terms found in the License Agreement file named 'COPYING' found in the
% FreeSurfer source code root directory, and duplicated here:
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
%
% General inquiries: freesurfer@nmr.mgh.harvard.edu
% Bug reports: analysis-bugs@nmr.mgh.harvard.edu
%


size = resx*resy;
fid = fopen(infile, 'r');
dummy = fread(fid, offset, type);
for i =1:(slice)
	A = fread(fid, size, type);
end
status = fclose(fid);

A = reshape(A,resx,resy);






