function err = fast_svcorslice(corslice,corslicefile)
%
% err = fast_svcorslice(corslice,corslicefile)
%


%
% fast_svcorslice.m
%
% Original Author: Doug Greve
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

if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE: err = fast_svcorslice(corslice,corslicefile)';
  qoe(msg);error(msg);
end

%corslice = fast_ldcorslice('tmp.cor');

% corslice = corslice'; %' Convert to column major 
corslice = uint8(corslice);

%%%% Open the corslicefile %%%%%
fid=fopen(corslicefile,'wb');
if(fid == -1)
  msg = sprintf('Could not open %s for writing.',corslicefile); 
  qoe(msg); error(msg);
end

%%% Write the file in corslicefile %%%
precision = 'uint8';
Nv = prod(size(corslice));
count = fwrite(fid,reshape1d(corslice),precision);
%count = fwrite(fid,reshape1d(corslice),'integer*1');
fclose(fid); 

if(count ~= Nv)
  fprintf(2,'ERROR: wrote %d/%d elements to %s\n',count,Nv,corslicefile);
  err = 1;
else 
  err = 0;
end

return;
