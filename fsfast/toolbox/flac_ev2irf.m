function Xirf = flac_ev2irf(flac,nthev)
% Xirf = flac_ev2irf(flac,nthev)
%
% $Id: flac_ev2irf.m,v 1.4 2006/10/20 03:39:58 greve Exp $

Xirf = [];
if(nargin ~= 2)
  fprintf('Xirf = flac_ev2irf(flac,nthev)\n');
  return;
end

ev = flac.ev(nthev);
if(~ev.ishrf)
  fprintf('ERROR: cannot get irf from non-HRF EV\n');
  return;
end

switch(ev.model)
 
 case {'fir'}
  Xirf = eye(ev.npsd);
 
 case {'spmhrf'}
  nderiv = ev.params(1);
  dpsd   = ev.params(2);
  t = dpsd*[0:ev.npsd-1]';
  Xirf = fast_spmhrf(t);
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
  t = dpsd*[0:ev.npsd-1]';
  Xirf = fmri_hemodyn(t,delay,tau,alpha);
  Xirf = Xirf/max(Xirf); % consistent with selxavg
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
  t = dpsd*[0:ev.npsd-1]' + flac.TR/2;
  Xirf = pdf_gamma(t,a,b);  
  Xirf = Xirf/sum(Xirf);
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









