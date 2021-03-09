function [cimg, Rrow, epipar] = tdr_epirecon(kepi,arg2,kpcn)
% [cimg Rrow epipar] = tdr_epirecon(kepi,measasc,<kpcn>)
% [cimg Rrow epipar] = tdr_epirecon(kepi,epipar,<kpcn>)
%
% Performs fouirer-based epi recon with ramp-sampling and 
% ghost correction.
%
% kepi is kspace data (nrow x nkcols x n1 x n2 x n3)
%  The interpretation of n1, n2, n3 is not important.
%  Applies readout reversals of even numbered rows.
% 
% EPI parameters are obtained in one of two ways:
%   measasc - a Siemens meas.asc 
%   epipar - a structure with the following fields:
%     tDelSamp  : delay from RF pulse to first readout adc (usec)
%     tDwell    : time to read one sample (usec)
%     tRampUp   : total ramp up time (usec)
%     tFlat     : total flat top time (usec)
%     tRampDown : total ramp down time (usec)
%     All times are in usec. Note: meas.asc has dwell time in nsec
%
% Returns an image with half the kspace cols.
% cimg is complex.
% 
% Reminder dont use "'" when transposing unless you 
% want complex conj.
%
% Example:
%  cd to dir with meas.asc and meas.out
%  mri_convert_mdh --srcdir . --binary --nopcn --binary --dumpadc adc.dat
%  fp = fopen('adc.dat','r');
%  ktmp = fread(fp,'float32');
%  fclose(fp);
%  kepi = ktmp(1:2:end) + i*ktmp(2:2:end); % real and imaginary
%  kepi = reshape(kepi,[nkcols nrows nchannels nslices ntp]);
%  kepi = permute(kepi,[2 1 3 4 5]);
%  [cimg Rrow epipar] = tdr_epirecon(kepi,'meas.asc');
%
%


%
% tdr_epirecon.m
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

cimg = [];
Rrow = [];
if(nargin < 2 | nargin > 3)
  fprintf('[cimg Rrow epipar] = tdr_epirecon(kepi,measasc,<kpcn>)\n');
  fprintf('[cimg Rrow epipar] = tdr_epirecon(kepi,epipar,<kpcn>)\n');
  return;
end

% The 2nd arg is either the path to a meas.asc file or an epipar struct
if(isstr(arg2))
  measasc = arg2;
  epipar.tDwell = tdr_measasc(measasc,'sRXSPEC.alDwellTime[0]'); % nsec
  if(isempty(epipar.tDwell)) return; end
  epipar.tDwell = epipar.tDwell/1000; % convert to usec
  epipar.tRampUp   = tdr_measasc(measasc,'alRegridRampupTime[0]');  % us
  epipar.tFlat     = tdr_measasc(measasc,'alRegridFlattopTime[0]'); % us
  epipar.tRampDown = tdr_measasc(measasc,'alRegridRampdownTime[0]'); % us
  epipar.tDelSamp  = tdr_measasc(measasc,'alRegridDelaySamplesTime[0]'); % us
  epipar.echospacing = tdr_measasc(measasc,'sFastImaging.lEchoSpacing'); % us
  epipar.nkcols    = tdr_measasc(measasc,'iNoOfFourierColumns'); 
  epipar.nkrows    = tdr_measasc(measasc,'iNoOfFourierLines'); 
else
  epipar = arg2;
end

if(~exist('kpcn','var')) kpcn = []; end

nrows   = size(kepi,1);
nkcols  = size(kepi,2);
% Higher dims, could be channels, slices, time points. 
% Does not really matter.
n1 = size(kepi,3); 
n2 = size(kepi,4);
n3 = size(kepi,5);

% Apply readout reversal to even numbered rows.
indkevenrows = [2:2:nrows];
kepi(indkevenrows,:,:,:) = flipdim(kepi(indkevenrows,:,:,:),2);

% Keep only the half of the image cols in the middle
c1 = nkcols/4 + 1;
c2 = c1 + nkcols/2 -1;
colkeep = [c1:c2];
ncolskeep = length(colkeep);

% Phase space trajectory for a single line (ramp-sampled)
kvec = kspacevector2(nkcols,epipar.tDwell,epipar.tRampUp,...
		     epipar.tFlat,epipar.tRampDown,epipar.tDelSamp,0);

Frow = fast_dftmtx(kvec);
Frow = Frow(:,colkeep);
Rrow = transpose(inv(Frow'*Frow)*Frow');

Fcol = fast_dftmtx(nrows);
Rcol = inv(Fcol);

% Recon the rows (note: tdr_deghost can sim)
simflag = 0;
kepi_rrecon = zeros(nrows,ncolskeep,n1,n2,n3);
for nth1 = 1:n1
  for nth2 = 1:n2
    for nth3 = 1:n3
      kimg = kepi(:,:,nth1,nth2,nth3);
      kepi_rrecon(:,:,nth1,nth2,nth3) = tdr_deghost(kimg,Rrow,0,simflag,kpcn);
    end
  end
end

% Recon the cols
kepi_rrecon = reshape(kepi_rrecon,[nrows ncolskeep*n1*n2*n3]);
cimg = Rcol*kepi_rrecon;
cimg = reshape(cimg, [nrows ncolskeep n1 n2 n3]);

return;

