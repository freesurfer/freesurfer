function [vavg, vstd, vjk] = fmri_jackknife(v,blflag)
% [vavg vstd] = fmri_jackknife(v,<blflag>)
%
% Computes the average and stddev of matrix v, where
% v is Nsamples x Nvar1 x Nvar2 ... 

if(nargin ~= 1 & nargin ~= 2)
  msg = 'USAGE: [vavg vstd] = fmri_jackknife(v,blflag)';
  qoe(msg);error(msg);
end

szv = size(v);
ndv = length(szv);
Nsamples = size(v,1);
Nvariables = prod(szv)/Nsamples;
v = reshape(v, [Nsamples Nvariables]);

for s = 1:Nsamples
  jk = find([1:Nsamples] ~= s);
  vjk(s,:) = mean(v(jk,:));
end

fudge = (Nsamples-1).^(1.5)/sqrt(Nsamples);
vavg = mean(vjk);
if(nargin == 2)
  vstd = fudge*sqrt(mean(vjk.^2));
else
  vstd = fudge*std(vjk);
end


if(ndv > 2)
  vavg = reshape(vavg, szv(2:ndv));
  vstd = reshape(vstd, szv(2:ndv));
end

return;
