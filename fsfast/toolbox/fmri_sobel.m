function [m, a] = fmri_sobel(img)
%
% [m, a] = fmri_sobel(img)
%

[dx dy] = gradient(img);
m = abs(dx + dy *sqrt(-1));


return;