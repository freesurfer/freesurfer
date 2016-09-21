function M = my_lta_read(fname)
% M = lta_read(fname)
% Modified by Eugenio so it reads lta files correctly when they do not have
% comments

if (strcmp(fname((length(fname)-2):length(fname)), 'xfm') | ...
	strcmp(fname((length(fname)-2):length(fname)), 'XFM'))
	M = xfm_read(fname) ;
	return ;
end


fid = fopen(fname) ;
if (fid < 0)
	error(sprintf('could not open file %s', fname));
end

tline = fgetl(fid) ;
while ~(length(tline) > 3  && strcmp(tline(1:4),'type')>0)
	tline = fgetl(fid) ;
end

tline = fgetl(fid) ;  % nxforms
tline = fgetl(fid) ;  % mean
tline = fgetl(fid) ;  % sigma
tline = fgetl(fid) ;  % dimensions

M = zeros(4,4) ;
for row=1:4
	tline = fgetl(fid) ;  % one row of matrix
	tmp = sscanf(tline, '%f');
	 M(row,:) = tmp';
end

fclose(fid) ;


