% fast_group_glm.m - see groupreg-sess
% $Id: fast_group_glm.m,v 1.1 2004/12/01 21:04:04 greve Exp $
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

tic;

ver = '$Id: fast_group_glm.m,v 1.1 2004/12/01 21:04:04 greve Exp $';
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

nvslice = nr*nc;

mristruct = fast_ldbhdr(instem);
if(isempty(mristruct))
  fprintf('ERROR: could not load bhdr for %s\n',instem);
  if(QuitOnError) exit; end
  return;
end

for slice = 1:ns
  fprintf('slice = %d (%g)\n',slice,toc);

  y = [];
  for n = 1:ninputs

    instem = deblank(InstemList(n,:));
    yn = fast_ldbslice(instem,slice-1);
    if(isempty(yn))
      fprintf('ERROR: loading %s\n',instem);
      if(QuitOnError) exit; end
      return;
    end

    yn = reshape(yn,[nvslice nf])';
    y = [y; yn];

    xflafile = deblank(FLAXMatList(n,:));
    Xflastruct = load(xflafile);
    Xfla = Xflastruct.Xfinal;
    NBetafla = size(Xfla,2);
    NBetatask = size(Cfla,2);
    Cfla0 = zeros(1,NBetafla);
    Cfla0(1:NBetatask) = Cfla;
    
  end

  [beta rvar] = fast_glmfit(y,X);
  [F, Fsig, ces] = fast_fratio(beta,X,rvar,C);

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
  tmp = -sign(ces).*log10(abs(Fsig));
  fast_svbslice(tmp,stem,slice-1,'',mristruct);

end

fmri_touch(okfile);
fprintf('matlab: fast_group_glm done (%g)\n',toc);

