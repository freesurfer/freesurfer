function [tsig, t, cesavg, cesstd, dof] = fast_igxavg2(n1,avg1,var1,n2,avg2,var2)
% [tsig, t, cesavg, cesstd, dof] = fast_igxavg2(n1,avg1,var1,n2,avg2,var2)
%
% Tests difference between means (small samples). avg1-avg2. The tsig
% is the signed log10 of the signficance and is the same size as the
% input. The t-test is single-sided.
%
% For formula, see:
%  Modern Elementary Statistics (8th), John E. Freund and Gary A. Simon,
%  1992, Prentice-Hall Englewood Cliffs, New Jersey. (Inside back cover).
%
% $Id: fast_igxavg2.m,v 1.1 2003/03/04 20:47:38 greve Exp $
%

tsig = [];
t = [];
cesavg = [];
cesstd = [];

if(nargin ~= 6)
  fprintf('[tsig t cesavg cesstd dof] = fast_igxavg2(n1,avg1,var1,n2,avg2,var2)\n');
  return;
end

dof = n1+n2-2;
cesavg = avg1-avg2;
cesstd = sqrt( ( (var1*(n1-1) + var2*(n2-1)) /dof ) * (1/n1 + 1/n2) );
indz = find(cesstd == 0);
cesstd(indz) = 1;

t = cesavg ./ cesstd;
t(indz) = 0;
tsig = tTest(dof,reshape1d(abs(t)));
tsig = tsig/2; % Single-sided

tsig = reshape(tsig, size(t));
tsig = sign(t) .* -log10(tsig);

return;
