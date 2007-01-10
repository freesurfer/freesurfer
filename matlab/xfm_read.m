function M = xfm_read(fname)
% M = xfm_read(fname)


%
% xfm_read.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:10 $
%    $Revision: 1.3 $
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


