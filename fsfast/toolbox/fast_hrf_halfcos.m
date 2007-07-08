function [hrf,params] = fast_hrf_halfcos(t,params)
% [hrf params] = fast_hrf_halfcos(t,<params>)
%
% Half-cosine HRF parameterization (used in FLOBS).
% 
% t - sample times.
% params = [h1 h2 h3 h4 f1 f2]
% h1, h2, h3, h4 - half period of segments
% f1 = max predip amp (at t=h1)
% f2 = max postdip amp (at t=h1+h2+h3)
% Peak is +1 at t=h1+h2 by construction. No scaling is performed.
%
% If params is not specified, then it is generated randomly
% according to the paper, almost. I had to multiply their periods
% by 1.5 to get something reasonable.
%
% Mark W. Woolrich, Timothy E.J. Behrens, and Stephen M. Smith.
% Constrained linear basis sets for HRF modelling using Variational
% Bayes NeuroImage 21 (2004) 1748 1761.
%
% $Id: fast_hrf_halfcos.m,v 1.2 2007/07/08 20:59:51 greve Exp $

hrf = [];
if(nargin < 1 & nargin > 2)
  fprintf('hrf = fast_hrf_halfcos(t,params)\n');
  return;
end

if(nargin == 2)
  % Parameters specified explicitly
  h1 = params(1);
  h2 = params(2);
  h3 = params(3);
  h4 = params(4);
  f1 = params(5);
  f2 = params(6);
else
  % Generate params randomly
  if(0)
    h1 = 2*rand;     % uniform(0,2), avg = 1
    h2 = 2 + 4*rand; % uniform(2,6), avg = 4
    h3 = 2 + 4*rand; % uniform(2,6), avg = 4 (too fast, 8 better)
    h4 = 2 + 6*rand; % uniform(2,8), avg = 6
    sc = 1.5; % Scale to make it look like the paper
    h1=sc*h1; h2=sc*h2; h3=sc*h3; h4=sc*h4;
    f1 = 0;
    f2 = .5*rand;    % uniform(0,.5), avg = .25, SPM=.89
  else
    h1 = 2*rand;         % uniform(0,2)
    h2 = 3 + (8-3)*rand; % uniform(3,8)
    h3 = 3 + (8-3)*rand; % uniform(3,6)
    h4 = 3 + (8-3)*rand; % uniform(3,8)
    f1 = 0;
    f2 = .3*rand;    % uniform(0,.3), avg = .3, SPM=.89
  end

  params = [h1 h2 h3 h4 f1 f2];
end

% Cosine Periods
P1 = 2*h1;
P2 = 2*h2;
P3 = 2*h3;
P4 = 2*h4;

% Segment onset times
t1  = 0;        
t2  = h1;       
t3  = h1+h2;    
t4  = h1+h2+h3; 
tend = h1+h2+h3+h4;

% Amplitudes
A1 = f1/2;
A2 = (1+f1)/2;
A3 = (1+f2)/2;
A4 = f2/2;

% Zero-crossings
z1 = -f1/2;
z2 = A2-f1;
z3 = A3-f2; 
z4 = -f2/2;

hrf = zeros(size(t));

ind1 = find(t > t1 & t <= t2);
hrf(ind1) = A1*cos(2*pi*(t(ind1)-t1)/P1) + z1;

ind2 = find(t > t2 & t <= t3);
hrf(ind2) = -A2*cos(2*pi*(t(ind2)-t2)/P2) + z2;

ind3 = find(t > t3 & t <= t4);
hrf(ind3) = A3*cos(2*pi*(t(ind3)-t3)/P3) + z3;

ind4 = find(t > t4 & t <= tend);
hrf(ind4) = -A4*cos(2*pi*(t(ind4)-t4)/P4) + z4;

return;







