function m = read_genesis_image(fname)
% m = read_ge(fname)

fid = fopen(fname, 'r', 'b');
if fid < 0
	disp(sprintf('read_genesis_image: error opening file %s', fname));
	return;
end

fseek(fid, 4, 'bof');
image_offset = fread(fid, 1, 'int');
rows = fread(fid, 1, 'int');
cols = fread(fid, 1, 'int');

fseek(fid, image_offset, 'bof');
m = fread(fid, [rows cols], 'short');

fclose(fid);

% eof
