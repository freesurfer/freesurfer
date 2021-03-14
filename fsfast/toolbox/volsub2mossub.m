function [rm, cm, tszmos] = volsub2mossub(rv, cv, sv, szvol, tszmos)
% [rm cm tszmos] = volsub2mossub(rv, cv, sv, szvol, tszmos)
%
% rv - row in the volume
% cv - column in the volume
% sv - slice in the volume
% szvol - size of the volume (Nrows, Ncols, Nslices, ...)
% tszmos - size (rows, cols) of the mosaic measured in tiles (optional)
%
% rm - row in the mosaic
% cm - column in the mosaic
%
% If tszmos is not specified, one will be computed.


%
% volsub2mossub.m
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

if(nargin ~= 4 & nargin ~= 5)
  msg = 'USAGE: [rm cm tszmos] = volsub2mossub(rv, cv, sv, szvol, <tszmos>)';
  error(msg);
end

if(length(rv) ~= length(cv) | ...
   length(rv) ~= length(sv) | ...
   length(cv) ~= length(sv) )
   msg = sprintf('rv (%d), cv (%d), and sv (%d) do not have the same length',...
                 length(rv),length(cv),length(sv));
   error(msg);
end

Nvr = szvol(1);
Nvc = szvol(2);
Nvs = szvol(3);

if(nargin == 4) tszmos = []; end
tszmos = defmossize(Nvs, tszmos);
Ntr = tszmos(1);
Ntc = tszmos(2);

rt = floor((sv-1)/Ntc) + 1;
ct = sv - (rt-1)*Ntc ;

rm = (rt-1)*Nvr + rv; 
cm = (ct-1)*Nvc + cv; 

%keyboard

return
