function m = read_siemens_image(fname)
% m = read_siemens(fname)

fid = fopen(fname, 'r', 'b');
if fid < 0
	disp(sprintf('read_siemens_image: error opening file %s', fname));
	return;
end

fseek(fid, 4994, 'bof');
rows = fread(fid, 1, 'short');
cols = fread(fid, 1, 'short');

fseek(fid, 6144, 'bof');
m = fread(fid, [rows cols], 'short');

fclose(fid);

% eof
