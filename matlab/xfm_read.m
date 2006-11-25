function M = xfm_read(fname)

% M = xfm_read(fname)

fid = fopen(fname) ;
if (fid < 0)
	error(sprintf('could not open file %s', fname));
end

tline = fgetl(fid) ;  
while ((length(tline) > 0) & (tline(1) == '%'))
	tline = fgetl(fid) ;
end

tok = strtok(tline);
while (strcmp(tok, 'Linear_Transform') ~= 1)
	tline = fgetl(fid) ;
	tok = strtok(tline);
end


M = zeros(4,4) ; M(4,4) = 1;
for row=1:3
	tline = fgetl(fid) ;  % one row of matrix
	tmp = sscanf(tline, '%f');
	 M(row,:) = tmp';
end

fclose(fid) ;


