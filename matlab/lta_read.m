function M = lta_read(fname)
% M = lta_read(fname)


%
% lta_read.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:12 $
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
while ((length(tline) > 0) & (tline(1) == '#'))
	tline = fgetl(fid) ;
end


tline = fgetl(fid) ;  % type
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


