function [s,fa] = ssbloch(tr,te,fa,t1,t2s,p)
% [s fa] = ssbloch(tr,te,fa,t1,t2s,<p>)
%
% Steady-state Bloch Equation
%
% tr = repetition time
% te = echo time
% fa = flip angle (radians)
% t1 = T1
% t2s = T2 star
% p = proton density (def = 1)
%
% Time units don't matter as long as they are consitent
%
% if fa=[], set to ernst angle:
%   fa = acos(exp(-TR/T1))
%
% From: Wansapura, et al, J MAG RESE IMG 9:531 538 (1999)
%  At 3T, 
%  Gray:  T1 = 1331ms, T2* = 42-52 ms
%  White: T1 =  832ms, T2* = 45-48 ms
%   

% $Id: ssbloch.m,v 1.1 2006/07/06 05:14:52 greve Exp $

s = [];
if(nargin < 5 | nargin > 6)
  fprintf('[s fa]= ssbloch(tr,te,fa,t1,t2s,<p>)\n');
  return;
end


if(~exist('p','var')) p = 1; end
if(isempty(fa)) fa = acos(exp(-tr./t1)); end

etrt1 = exp(-tr./t1);
s = p .* sin(fa) .* (1-etrt1) .* exp(-te./t2s) ./ (1-cos(fa).*etrt1);

return;






