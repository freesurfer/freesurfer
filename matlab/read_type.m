function m = read_type(fname, offset, size, type)
% m = read_type('fname', offset, size, 'type')


%
% read_type.m
%
% Original Author: Bruce Fischl
%
% Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
%
% Terms and conditions for use, reproduction, distribution and contribution
% are found in the 'FreeSurfer Software License Agreement' contained
% in the file 'LICENSE' found in the FreeSurfer distribution, and here:
%
% https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
%
% Reporting: freesurfer@nmr.mgh.harvard.edu
%

fid = fopen(fname, 'r', 'b');
if fid < 0
	disp(sprintf('read_type: error opening file %s', fname));
	return;
end

fseek(fid, offset, 'bof');
m = fread(fid, size, type);

fclose(fid);

% eof
