function [Xirf tirf] = flac_ev2irf(ev,TR,RefEventDur)
% [Xirf tirf] = flac_ev2irf(ev,TR,RefEventDur)
%
% ev - explanatory variable with params defined
% TR - TR in seconds
% RefEventDur - Duration of reference event in seconds. This 
% sets the scaling factor. One can expect an event of this 
% duration witha peak of 1% signal change to have a beta=1.
%
% $Id: flac_ev2irf.m,v 1.8 2010/03/25 18:50:00 greve Exp $

%
% flac_ev2irf.m
%
% Original Author: Doug Greve
% CVS Revision Info:
%    $Author: greve $
%    $Date: 2010/03/25 18:50:00 $
%    $Revision: 1.8 $
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

Xirf = [];
if(nargin ~= 3)
  fprintf('[Xirf tirf] = flac_ev2irf(ev,TR,RefEventDur)\n');
  return;
end

if(~ev.ishrf)
  fprintf('ERROR: cannot get irf from non-HRF EV\n');
  return;
end

switch(ev.model)
 
 case {'fir'}
  Xirf = eye(ev.npsd);
  tirf = [];
 
 case {'spmhrf'}
  nderiv = ev.params(1);
  dpsd   = ev.params(2);
  tirf = dpsd*[0:ev.npsd-1]';
  Xirf = fast_spmhrf_sampled(tirf);
  % Scale for reference event duration
  Nref = round(RefEventDur/dpsd);
  a = ones(Nref,1);
  c = conv(Xirf,a);
  scalef = max(c);
  Xirf = Xirf/scalef;
  dhspmhrf = Xirf;
  for n = 1:nderiv
    % Divide by TER for gradient.
    % Multiply by 2.6 to bring 1st deriv to amp of 1
    dhspmhrf = 2.6*gradient(dhspmhrf)/dpsd;
    Xirf = [Xirf dhspmhrf];
  end  

 case {'gamma'}
  delay  = ev.params(1);
  tau    = ev.params(2);
  alpha  = ev.params(3);
  nderiv = ev.params(4);
  dpsd   = ev.params(5);
  tirf = dpsd*[0:ev.npsd-1]';
  Xirf = fmri_hemodyn(tirf,delay,tau,alpha);
  % Scale for reference event duration
  Nref = round(RefEventDur/dpsd);
  a = ones(Nref,1);
  c = conv(Xirf,a);
  scalef = max(c);
  Xirf = Xirf/scalef;
  %Xirf = Xirf/max(Xirf); % this is kinda what selxavg does
  dh_hrf = Xirf;
  for n = 1:nderiv
    % Divide by TER for gradient.
    % Multiply by 2.6 to bring 1st deriv to amp of 1
    dh_hrf = 2.6*gradient(dh_hrf)/dpsd;
    Xirf = [Xirf dh_hrf];
  end  

 case {'fslgamma'}
  phase   = ev.params(1); % Not used yet
  sigma   = ev.params(2);
  meanlag = ev.params(3);
  a = (meanlag/sigma)^2;
  b = meanlag/(sigma^2);
  nderiv  = ev.params(4);
  dpsd    = ev.params(5);
  tirf = dpsd*[0:ev.npsd-1]' + TR/2;
  Xirf = pdf_gamma(tirf,a,b);  
  % Scale for reference event duration
  Nref = round(RefEventDur/dpsd);
  a = ones(Nref,1);
  c = conv(Xirf,a);
  scalef = max(c);
  Xirf = Xirf/scalef;
  %Xirf = Xirf/sum(Xirf);
  dh_hrf = Xirf;
  for n = 1:nderiv
    % Divide by TER for gradient.
    % Multiply by 2.6 to bring 1st deriv to amp of 1
    dh_hrf = 2.6*gradient(dh_hrf)/dpsd;
    Xirf = [Xirf dh_hrf];
  end  

 otherwise 
  fprintf('ERROR: model %s unrecognized\n',ev.model);
  
end

return;









