function M = lta_read(fname)
% M = lta_read(fname)


%
% lta_read.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: nicks $
%    $Date: 2007/01/10 22:55:09 $
%    $Revision: 1.2 $
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


