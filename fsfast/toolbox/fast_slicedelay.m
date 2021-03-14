function delay = fast_slicedelay(TR,nSlices,Slice,AcqOrder)
% delay = fast_slicedelay(TR,nSlices,Slice,AcqOrder)
%
% Computes the amount of time after the start of a TR that a slice
% is aquired given the TR, the number of slices, and the order of
% acquisition (linear or interleaved).
%


%
% fast_slicedelay.m
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

if(nargin ~= 4)
     msg = 'USAGE: delay = fast_slicedelay(TR,nSlices,Slice,AcqOrder)';
  qoe(msg);error(msg);
end

if( ~strcmpi(AcqOrder,'linear') & ~strcmpi(AcqOrder,'interleaved'))
  msg = sprintf('AcqOrder = %s, must be either linear or interleaved',AcqOrder);
  qoe(msg);error(msg);
end

dt = TR/nSlices;

if( strcmpi(AcqOrder,'linear') )
  SliceOrder = [0:nSlices-1];
end

if( strcmpi(AcqOrder,'interleaved') )
  SliceOrder = [[0:2:nSlices-1] [1:2:nSlices] ];
end

nthAcq = find(SliceOrder == Slice) - 1;
delay = dt*nthAcq;

return; 
