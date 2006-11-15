% fast_selxavg3.m
% $Id: fast_selxavg3.m,v 1.1 2006/11/15 00:31:49 greve Exp $

% analysis
% sessdir

monly = 1;
analysis = 'edp2';
sess  = 'tl20000621';

flac0 = fast_ldanaflac(analysis);
if(isempty(flac0))
  if(~monly) quit; end
  return; 
end

flac0.sess = sess;
flac0.nthrun = 1;
flac0 = flac_customize(flac0);
if(isempty(flac0)) 
  if(~monly) quit; end
  return; 
end

nruns = size(flac0.runlist,1);
fprintf('nruns = %d\n',nruns);

fprintf('Creating Design Matrix\n');
Xt = [];
Xn = [];
tic;
for nthrun = 1:nruns
  flac = flac0;
  flac.nthrun = nthrun;
  flac = flac_customize(flac);
  if(isempty(flac)) 
    if(~monly) quit;  end
    return; 
  end

  indtask = flac_taskregind(flac);
  Xtr = flac.X(:,indtask);

  indnuis = flac_nuisregind(flac);
  Xnr = flac.X(:,indnuis);

  Za = zeros(size(Xn,1),  size(Xnr,2));
  Zb = zeros(size(Xnr,1), size(Xn,2));
  
  Xt = [Xt; Xtr];
  Xn = [Xn Za; Zb Xnr];
  
end

X = [Xt Xn];






