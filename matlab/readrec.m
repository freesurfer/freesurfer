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

size = resx*resy;
fid = fopen(infile, 'r');
dummy = fread(fid, offset, type);
for i =1:(slice)
	A = fread(fid, size, type);
end
status = fclose(fid);

A = reshape(A,resx,resy);






