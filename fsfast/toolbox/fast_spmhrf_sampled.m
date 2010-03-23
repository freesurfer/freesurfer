function [hrf, tsamples] = fast_spmhrf_sampled(tsamples,nderiv)
% [hrf tsamples] = fast_spmhrf_sampled(tsamples,<nderiv>)
% 
% Samples spmhrf at the given times.  fast_spmhrf will try
% to do some normalization, which changes the shape/amplitude when
% sampled over different times.  This works by sampling the HRF at a
% very fine resolution (1ms), then finding the values at tsamples in
% the oversampled HRF. The HRF is scaled such that the
% *oversampled* HRF is 1.0. The max of the returned hrf may or may
% not be 1.0 depending upon how close one of the tsamples falls
% to where the max is. 
%
% If nderiv is > 0, then the gradient of the hires HRF is computed
% nderiv times then sampled at tsample. The 1ms sample time is
% factored in when compuding the gradient: hrf0 = gradient(hrf0,dt);
%
% tsamples can be a matrix.
%
% TR = 2;
% nslices = 30;
% dtslice = TR/nslices;
% tslice = dtslice*[0:nslices-1];
% tfir = [0:TR:32]';
% nfir = length(tfir);
% tsamples = repmat(tfir,[1 nslices]) - repmat(tslice,[nfir 1]);
% hrf = fast_spmhrf_sampled(tsamples);
% hrf = hrf/sum(hrf(:,1));
%
% 
% $Id: fast_spmhrf_sampled.m,v 1.2 2010/03/23 17:00:21 greve Exp $

hrf = [];
if(nargin < 1 | nargin > 2)
  fprintf('[hrf tsamples] = fast_spmhrf_sampled(tsamples,<nderiv>)\n');
  return;
end

if(nargin == 1) nderiv = []; end
if(isempty(nderiv)) nderive = 0; end

% Sample HRF at a very fine resolution
dt = .001;

tsamplesmax = dt*round(max(tsamples(:))/dt);
tsamplesmin = dt*round(min(tsamples(:))/dt);
twindow = [tsamplesmin:dt:tsamplesmax]';

hrf0 = fast_spmhrf(twindow);
hrf0 = hrf0/max(hrf0);
if(nderiv > 0)
  for n = 1:nderiv
    hrf0 = gradient(hrf0,dt);
  end
end

isamples = round((tsamples-tsamplesmin)/dt) + 1;
hrf = hrf0(isamples);

tsamples = dt*(isamples-1) + tsamplesmin;

return;





