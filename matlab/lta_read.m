function M = lta_read(fname)
% M = lta_read(fname)


%
% lta_read.m
%
% Original Author: Bruce Fischl
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2014/07/10 04:11:44 $
%    $Revision: 1.5 $
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

n=0;
while(1)
  n = n + 1;
  tline = fgetl(fid);
  if(tline == -1)
    fprintf('lta_read format error %s\n',fname);
    M = [];
    return;
  end
  tag = sscanf(tline,'%s',1);
  if(strcmp(tag,'type')) break; end
end

tline = fgetl(fid);   % nxforms
tline = fgetl(fid);   % mean
tline = fgetl(fid);   % sigma
tline = fgetl(fid);   % dimensions

M = zeros(4,4) ;
for row=1:4
  tline = fgetl(fid);   % one row of matrix
  tmp = sscanf(tline, '%f');
  M(row,:) = tmp';
end

fclose(fid) ;


