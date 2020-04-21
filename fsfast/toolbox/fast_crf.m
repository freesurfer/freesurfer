function crf = fast_crf(t)
% crf = fast_crf(t)
%
% Cardiac Response Function from Chang, et al, NI 44 (2009),
% 857-869, Equation 5 and Figure 6.
%

if(nargin ~= 1)
  fprintf('crf = fast_crf(t)\n');
  return;
end

e = -0.5*((t-12).^2)/9;

crf = 0.6 * (t.^2.7) .* exp(-t/1.6) - 16*exp(e)/sqrt(2*pi*9);

return;











