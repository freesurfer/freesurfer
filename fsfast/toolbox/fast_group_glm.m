% fast_group_glm.m - see groupreg-sess
%
%
% InstemList = splitstring('$InstemList');
% FLAXMatList = splitstring('$FLAXMatList');
% FLAConMat = '$FLAConMat';
% ninputs = size(InstemList,1);
% outdir = '$outdir';
% okfile = '$okfile';
% xmatfile = '$xmat';
% gconmatfile = '$gconmat';
% QuitOnError = ~[$monly];
% hemicode = 'lh'; hemicode = '';
% nthframe = 1;


%
% fast_group_glm.m
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

tic;

ver = 'fast_group_glm.m @FS_VERSION@';
fprintf('%s\n',ver);

Cflastruct = load(FLAConMat);
if(isempty(Cflastruct))
  fprintf('ERROR: %s\n',FLAConMat);
  return;
end
Cfla = Cflastruct.ContrastMtx_0;

% X = load(xmatfile,'-ascii');
[X, rowid, colid] = fast_ldtable(xmatfile);
if(isempty(X))
  fprintf('ERROR: loading %s\n',xmatfile);
  return;
end
[nxf nbeta] = size(X);

XtX = X'*X;
c = sqrt(cond(XtX));
fprintf('Design matrix condition %g\n',c);
if(c > 1000)
  fprintf('ERROR: design matrix is ill-conditioned.\n');
  fprintf('  Run groupreg-sess -help for more info.\n');
  if(QuitOnError) exit; end
  return;
end

C = load(gconmatfile,'-ascii');
if(isempty(C))
  fprintf('ERROR: loading %s\n',gconmatfile);
  if(QuitOnError) exit; end
  return;
end
ncbeta = size(C,2);

if(nbeta ~= ncbeta)
  fprintf('ERROR: number of columns in X (%d) does not equal\n');
  fprintf('  the number of columns in C (%d)\n',nbeta,ncbeta);
  if(QuitOnError) exit; end
  return;
end

if(nxf ~= ninputs)
  fprintf('ERROR: number of rows in X (%d) does not equal\n');
  fprintf('  the number of input (%d)\n',nxf,ninputs);
  if(QuitOnError) exit; end
  return;
end

instem = deblank(InstemList(1,:));
[ns nr nc nf] = fmri_bvoldim(instem);
if(isempty(ns))
  fprintf('ERROR: loading %s\n',instem);
  return;
end
if(nf > 1)
  fprintf('ERROR: first-level contrast has more than one plane\n');
  if(QuitOnError) exit; end
  return;
end

if(nthframe > nf)
  fprintf('ERROR: frame=%d > nframes=%d\n',nthframe,nf);
  return;
end

nvslice = nr*nc;

mristruct = fast_ldbhdr(instem);
if(isempty(mristruct))
  fprintf('ERROR: could not load bhdr for %s\n',instem);
  if(QuitOnError) exit; end
  return;
end

if(synth)
  fprintf('Synthsizing input data\n');
  SynthSeed = sum(100*clock);
  fprintf('SynthSeed = %10d\n',SynthSeed);
end

for slice = 1:ns
  fprintf('slice = %d (%g)\n',slice,toc);

  y = [];
  yvar = [];
  for n = 1:ninputs

    instem = deblank(InstemList(n,:));
    yn = fast_ldbslice(instem,slice-1);
    if(isempty(yn))
      fprintf('ERROR: loading %s\n',instem);
      if(QuitOnError) exit; end
      return;
    end
    yn = yn(:,:,nthframe);
    yn = reshape(yn,[nvslice 1])';

    if(synth) yn = randn(size(yn)); end
    
    if(WLS)
      invarstem = deblank(InVarStemList(n,:));      
      ynvar = fast_ldbslice(invarstem,slice-1);
      if(isempty(ynvar))
	fprintf('ERROR: loading %s\n',invarstem);
	if(QuitOnError) exit; end
	return;
      end
      ynvar = ynvar(:)';
      yvar = [yvar; ynvar];
      if(synth) yn = yn .* sqrt(ynvar); end
    end
    
    y = [y; yn];

    xflafile = deblank(FLAXMatList(n,:));
    Xflastruct = load(xflafile);
    Xfla = Xflastruct.Xfinal;
    NBetafla = size(Xfla,2);
    NBetatask = size(Cfla,2);
    Cfla0 = zeros(1,NBetafla);
    Cfla0(1:NBetatask) = Cfla;
    
  end % Loop over inputs
  
  inddata = find(all(y));
  ndata = length(inddata);
  indnodata = find(any(y==0));
  fprintf('  Found %d voxels with data for all inputs\n',ndata);
  if(WLS)
    yvarsum = sum(yvar);
    nbeta = size(X,2);
    beta = zeros(nbeta,nvslice);
    rvar = zeros(1,nvslice);
    F    = zeros(1,nvslice);
    Fsig = ones(1,nvslice);
    ces  = zeros(1,nvslice);
    % weight is 1/std
    w = zeros(ninputs,nvslice);
    w(:,inddata) = 1./sqrt(yvar(:,inddata));
    % Rescale weight so that the sum=1 at each voxel
    wsum = sum(w);
    w(:,inddata) = w(:,inddata) ./ repmat(wsum(inddata),[ninputs 1]);
    X0 = X;
    fprintf('  Staring WLS loop over %d voxels (%g)\n',length(inddata),toc);
    for nthind = 1:ndata
      indv = inddata(nthind); % index into the volume
      wv = w(:,indv);
      yv = wv.*y(:,indv);
      Xv = X .* repmat(wv,[1 nbeta]);
      [betav rvarv] = fast_glmfit(yv,Xv);
      [Fv, Fsigv, cesv] = fast_fratio(betav,Xv,rvarv,C);
      beta(:,indv) = betav;
      rvar(indv) = rvarv;
      F(indv) = Fv;
      Fsig(indv) = Fsigv;
      ces(indv) = cesv;
    end
  else
    [beta rvar] = fast_glmfit(y,X);
    [F, Fsig, ces] = fast_fratio(beta,X,rvar,C);
    beta(:,indnodata) = 0;
    rvar(indnodata) = 0;
    F(indnodata) = 0;
    Fsig(indnodata) = 1;
    ces(indnodata) = 1;
  end

  beta = reshape(beta', [nr nc nbeta]);
  rvar = reshape(rvar', [nr nc 1]);
  F    = reshape(F',    [nr nc 1]);
  Fsig = reshape(Fsig', [nr nc 1]);
  ces  = reshape(ces',  [nr nc 1]);

  stem = sprintf('%s/beta%s',outdir,hemicode);
  fast_svbslice(beta,stem,slice-1,'',mristruct);

  stem = sprintf('%s/beta-var%s',outdir,hemicode);
  fast_svbslice(rvar,stem,slice-1,'',mristruct);

  stem = sprintf('%s/f%s',outdir,hemicode);
  fast_svbslice(F,stem,slice-1,'',mristruct);

  stem = sprintf('%s/ces%s',outdir,hemicode);
  fast_svbslice(ces,stem,slice-1,'',mristruct);

  stem = sprintf('%s/sig%s',outdir,hemicode);
  tmp = zeros(size(Fsig));
  indok = find(Fsig ~= 0);
  tmp(indok) = -sign(ces(indok)).*log10(abs(Fsig(indok)));
  fast_svbslice(tmp,stem,slice-1,'',mristruct);

end

fmri_touch(okfile);
fprintf('matlab: fast_group_glm done (%g)\n',toc);

