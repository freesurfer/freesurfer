function m = read_type(fname, offset, size, type)
% m = read_type('fname', offset, size, 'type')

fid = fopen(fname, 'r', 'b');
if fid < 0
	disp(sprintf('read_type: error opening file %s', fname));
	return;
end

fseek(fid, offset, 'bof');
m = fread(fid, size, type);

fclose(fid);

% eof
