function M = xfm_read(fname)
% M = xfm_read(fname)


%
% xfm_read.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2011/03/02 00:04:13 $
%    $Revision: 1.4 $
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


