function Xirf = flac_ev2irf(flac,nthev)
% Xirf = flac_ev2irf(flac,nthev)
%
% $Id: flac_ev2irf.m,v 1.1 2004/10/17 18:34:48 greve Exp $

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

 otherwise 
  fprintf('ERROR: model %s unrecognized\n',ev.model);
  
end

return;









