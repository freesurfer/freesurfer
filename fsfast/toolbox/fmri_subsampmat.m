function [Mss, Rss] = fmri_subsampmat(TR,Ntp,TS)
% [Mss Rss] = fmri_subsampmat(TR,Ntp,TS)
%
% Creates a subsampling matrix


if(nargin ~= 3)
  msg = 'USAGE: [Mss Rss] = fmri_subsampmat(TR,Ntp,TS)';
  qoe(msg);error(msg);
end

Rss = floor(TR/TS);
if(Rss < 1)
  msg = sprintf('Subsampling factor must be >= 1');
  qoe(msg);error(msg);
end

if(Rss == 1)
  Mss = eye(Ntp);
  return;
end

Nsp = Rss*Ntp;

Mss = zeros(Ntp,Nsp);
w   = [1:Rss:Nsp];
z   = [1:Ntp];
ind = sub2ind([Ntp Nsp], z,w);
Mss(ind) = 1;


return;
