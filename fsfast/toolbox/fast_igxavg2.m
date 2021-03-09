function [tsig, t, cesavg, cesstd, dof] = fast_igxavg2(n1,avg1,var1,n2,avg2,var2)
% [tsig, t, cesavg, cesstd, dof] = fast_igxavg2(n1,avg1,var1,n2,avg2,var2)
%
% Tests difference between means (small samples). avg1-avg2. The tsig
% is the signed log10 of the signficance and is the same size as the
% input. The t-test is single-sided.
%
% avg1 - average  of sample 1
% var1 - variance of sample 1 (not the variance of avg1)
% avg2 - average  of sample 2
% var2 - variance of sample 2 (not the variance of avg2)
%
% For formula, see:
%  Modern Elementary Statistics (8th), John E. Freund and Gary A. Simon,
%  1992, Prentice-Hall Englewood Cliffs, New Jersey. (Inside back
%  cover and page 324).
%
%
%


%
% fast_igxavg2.m
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
