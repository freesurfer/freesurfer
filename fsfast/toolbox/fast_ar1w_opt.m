function [alpha, rho, mse, niters] = fast_ar1w_opt(Mr,R,tol,ncycles,alpharho0)
% [alpha, rho, mse, niters] = fast_ar1w_opt(Mr,R,tol,ncycles,alpharho0)


%
% fast_ar1w_opt.m
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

alpha = [];
rho = [];
e = [];

if(nargin ~= 4 & nargin ~= 5)
  fprintf('[alpha, rho, e] = fast_ar1w_opt(Mr,R,tol,ncycles,<alpharho0>)\n');
  return;
end

if(0)
niters = 0;
mse = fast_ar1w_mse(.1,.01,Mr,R);
for alpha0 = .1:.1:.8
  for rho0 = .01:.01:.3
    msetry = fast_ar1w_mse(alpha0,rho0,Mr,R);
    if(mse > msetry) 
       mse = msetry;
       alpha = alpha0;
       rho = rho0;
    end
    niters = niters + 1;
  end
end
return;
end

if(nargin == 5)
  alpha = alpharho0(1);
  rho = alpharho0(2);
else
  alpha = 0.5;
  rho = 0.5;
end

fprintf('Starting at: alpha = %g, rho = %g\n',alpha, rho);

tic;
niters = 0;
for cycle = 1:ncycles
  alpha0 = alpha;
  rho0 = rho;
  for minwhich = 1:2
    if(minwhich == 1)
       [alpha aminerr aniters] = goldenmin(Mr,R,minwhich,rho,tol);
    else
       [rho rminerr rniters] = goldenmin(Mr,R,minwhich,alpha,tol);
    end
  end
  niters = niters + aniters + rniters;
  fprintf('cycle = %d %6.4f %6.4f %g %g %g %g \n',...
        cycle,alpha,rho,aminerr,rminerr,fast_ar1w_mse(alpha,rho,Mr,R),toc);
  if(alpha0-alpha == 0 & rho0-rho == 0)
    fprintf('INFO: No change in parameters\n');
    break;
  end
  if(rho == 0)
    fprintf('INFO: No rho = 0, breaking\n');
    break;
  end

end

mse = fast_ar1w_mse(alpha,rho,Mr,R);

return;

%--------------------------------------------------------%
function [v, minerr, niters] = goldenmin(Mr,R,minwhich,other,tol)

a = 0;
c = 1-tol;
b = a + (1-.618)*(c-a);

if(minwhich == 1)
  ya = fast_ar1w_mse(a,other,Mr,R);
  yb = fast_ar1w_mse(b,other,Mr,R);
  yc = fast_ar1w_mse(c,other,Mr,R);
else
  ya = fast_ar1w_mse(other,a,Mr,R);
  yb = fast_ar1w_mse(other,b,Mr,R);
  yc = fast_ar1w_mse(other,c,Mr,R);
end

niters = 0;
while(abs(a-c) > tol)
  dba = b-a;
  dcb = c-b;
  if(dba > dcb) x = a + .618*dba;
  else          x = b + .618*dcb;
  end
  if(minwhich == 1) yx = fast_ar1w_mse(x,other,Mr,R);
  else              yx = fast_ar1w_mse(other,x,Mr,R);
  end
  %fprintf('%3d a=%5.4f, b=%5.4f, c=%5.4f, x=%5.4f   %g\n',...
  %       niters,a,b,c,x,abs(a-c));
  %fprintf('%3d ya=%g, yb=%g, yc=%g, yx=%g \n',niters,ya,yb,yc,yx);
  if(yb < yx & b < x)     c=x; yc=yx;
  elseif(yb < yx & b > x) a=x; ya=yx;
  elseif(yb > yx & b < x) a=b; ya=yb; b=x; yb=yx; 
  elseif(yb > yx & b > x) c=b; yc=yb; b=x; yb=yx;
  end
  niters = niters + 1;
end

y = [ya yb yc];
abc = [a b c];

[minerr i] = min(y);
v = abc(i);
return;
%--------------------------------------------------------%




