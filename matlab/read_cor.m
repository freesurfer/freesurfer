function a = read_cor(fname)
% a = read_cor(fname)
a = zeros(256, 256, 256);

for i=1:256
%	disp(i);
	fname1 = sprintf('%s/COR-%03d', fname, i);
	fid = fopen(fname1, 'r');
	a(:, :, i) = fread(fid, [256 256], 'uchar');
	fclose(fid);
end

% eof
