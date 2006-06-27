function err = fast_xmat2xtxt(matfile,txtfile,omnimtxfile)
% err = fast_xmat2xtxt(matfile,txtfile,<omnimtxfile>)
%
% Converts an fsfast/selxavg design matrix into a text file 
% that can be read by mri_glmfit. Optionally also creates
% an omnibus contrast matrix.
% 
%   matfile - X.mat as created by selxavg
%   txtfile - text file with design matrix
%   omnimtxfile - contrast to test omnibus
% 
% Example:
% cd bold 
% In matlab
%   fast_xmat2xtxt('analysis/X.mat','x.txt','omnibus.mat')
% In unix: concat all the runs
%   mri_concat [0-9][0-9][0-9]/fmcsm5.bhdr --o cc.mgh 
% Run mri_glmfit (this may take 5-10 min)
%   mri_glmfit --glmdir glmfit --y cc.mgh --X x.txt --C omnibus.mat
% Look at omnibus results. Should be same as omnibus/fsig.bhdr in
% fsfast, but the fsfast version will be unsigned.
%   tkmedit subjectname orig.mgz -overlay-reg register.dat \
%      -overlay glmfit/omnibus/sig.mgh
%
% $Id: fast_xmat2xtxt.m,v 1.1 2006/06/27 19:52:59 greve Exp $

err = 1;

if(nargin < 2 & nargin > 3)
  fprintf('err = fast_xmat2xtxt(matfile,txtfile,<omnimtxfile>)\n');
  return;
end

tmp = load(matfile);
if(~exist('tmp','var'))
  fprintf('ERROR: loading %s\n',matfile);
  return;
end
X = tmp.Xfinal;
nreg = size(X,2);

fp = fopen(txtfile,'w');
if(fp == -1)
  fprintf('ERROR: opening %s for writing\n',txtfile);
  return;
end

fmt = [repmat('%g ',[1 nreg]) '\n'];
fprintf(fp,fmt,X');
fclose(fp);

if(nargin == 3)
  fp = fopen(omnimtxfile,'w');
  if(fp == -1)
    fprintf('ERROR: opening %s for writing\n',omnimtxfile);
    return;
  end
  ntaskreg = tmp.Navgs_per_cond*tmp.Nnnc;
  C = eye(ntaskreg,nreg);
  fprintf(fp,fmt,C');
  fclose(fp);
end


err = 0;

return
