function flacnew = flac_customize(flac)
%  flacnew = flac_customize(flac)
%
% $Id: flac_customize.m,v 1.1 2004/10/16 05:19:36 greve Exp $

flacnew = [];

if(nargin ~= 1)
 fprintf('[X, taskind] = flac_customize(flac)\n');
 return;
end

sess = flac.sess;
run  = flac.run;
runpath = sprintf('%s/%s/%03d',sess,flac.fsd,run);
fstem = sprintf('%s/%s',runpath,flac.funcstem);

[nslices nrows ncols ntp] = fmri_bvoldim(fstem);
if(nslices == 0) 
  fprintf('ERROR: attempting to read %s\n',fstem);
  return; 
end
flacnew = flac;
flacnew.ntp = ntp;

nev = length(flac.ev);
for nthev = 1:nev
  ev = flac.ev(nthev);
  
  % HRF Regressors
  if(ev.ishrf)  
    stfpath = sprintf('%s/%s',runpath,ev.stf);
    st = fast_ldstf(stfpath);
    if(isempty(st)) flacnew = []; return; end
    flacnew.ev(nthev).st = st;
    ev.Xfir = fast_st2fir(st,ntp,flac.TR,ev.psdwin,1);
  end  
  
  % Non-parametric regressors
  if(strcmp(ev.model,'nonpar'))
    nonparpath = sprintf('%s/%s',runpath,ev.nonparname);
    X = fast_ldbslice(nonparpath);
    if(isempty(X))
      fprintf('ERROR: loading %s\n',nonparpath);
      flacnew = [];
      return;
    end
    if(size(X,2)~=1) X = squeeze(X)';
    else             X = squeeze(X);
    end
    if(size(X,1) ~= ntp)
      fprintf('ERROR: nonpar time point mismatch %s\n',nonparpath);
      flacnew = [];
      return;
    end
    if(size(X,2) < ev.params(1))
      fprintf('ERROR: not enough columns %s\n',nonparpath);
      flacnew = [];
      return;
    end
    % Demean
    X = X(:,1:ev.params(1));
    Xmn = mean(X,1);
    X = X - repmat(Xmn,[ntp 1]);
    flacnew.ev(nthev).X = X;
  end
  
end

return;














