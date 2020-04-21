function rrf = fast_rrf(t)
% rrf = fast_rrf(t)
%
% Respiratory Response Function from Birn, et al, NI 40 (2008),
% 644-654, equation 3.
%

if(nargin ~= 1)
  fprintf('rrf = fast_rrf(t)\n');
  return;
end

rrf = .6 * (t.^2.1) .* exp(-t/1.6) - .0023 * (t.^3.54) .* exp(-t/4.25);

return;











