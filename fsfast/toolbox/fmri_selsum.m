function ss = fmri_selsum(X,y,Cn)
%
% ss = fmri_selsum(X,y,<Cn>)
%
%
% $Id: fmri_selsum.m,v 1.1 2003/03/04 20:47:40 greve Exp $
%

if(nargin ~= 2 & nargin ~= 3)
  msg = 'USAGE: ss = fmri_selsum(X,y,<Cn>)';
  qoe(msg);
  error(msg);
end

ss = 0;

[nTP nTotEst nRuns] = size(X);
[nRows nCols nTP nRuns2] = size(y);

if(nRuns ~= nRuns2)
  msg = 'X and y dimensions are inconsistent';
  qoe(msg);
  error(msg);
end

nV = nRows*nCols;

y = reshape(y, [nV nTP nRuns]);
y = permute(y, [2 1 3]);

if(nargin == 2) % Assume white noise, ie Cn = I
  for r = 1:nRuns,
     ss = ss + X(:,:,r)' * y(:,:,r); %'
  end

else % Use Cn as noise covariance matrix

  if(size(Cn,1) ~= nTP | size(Cn,2) ~= nTP | ...
     (size(Cn,3) ~=1 & size(Cn,3) ~= nRuns))
    msg = 'X and Cn dimensions are inconsistent';
    qoe(msg);
    error(msg);
  end

  if(size(Cn,3) == 1) % only one matrix specified, use for all runs
    iCn = inv(Cn);
    for r = 1:nRuns,
       ss = ss + X(:,:,r)' * iCn * y(:,:,r);
    end
  else  % use different Cn for each run
    for r = 1:nRuns,
       ss = ss + X(:,:,r)' * inv(Cn(:,:,r)) * y(:,:,r);
    end  
  end

end

nTotEst = size(X,2);
ss = permute(ss, [2 1 3]);
ss = reshape(ss, [nRows nCols nTotEst]);

return
