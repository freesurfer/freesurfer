function [err irepistruct] = irepierr(params,irepistruct)
% [err irepiparstruct] = irepierr(params,irepistruct)
% if(s.paramset == 1) T1 = params(1);
%

s = irepistruct; % copy into s for easy handling

switch s.parset
  case 1
   s.T1 = params(1);
  case 2
   s.T1 = params(1);
   s.sigma = params(2);
end

irepistruct = s;
if(s.T1 < 0  | s.T1 > 10000) err = 10^10; return; end
if(s.eff < 0 | s.eff > 1.5) err = 10^10; return; end
if(s.sigma < 0) err = 10^10; return; end

s = irepifit(s);
%err = s.rstd;
err = mean(abs(s.res));

irepistruct = s;
return;









