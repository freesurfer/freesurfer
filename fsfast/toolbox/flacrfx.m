% flacrfx.m
% $Id: flacrfx.m,v 1.1 2005/01/23 17:13:12 greve Exp $
tic;

%sess = '/homes/4/greve/links/sg1/skb/skb-101.1';
%flacdir = '/homes/4/greve/links/sg1/skb/flac/';
%flacfile = sprintf('%s/fn-swf-nvr.flac',flacdir);

%sess = '/homes/4/greve/links/sg1/dng081403/dng-siemens';
%flacdir = '/homes/4/greve/links/sg1/dng081403/flac';
%flacfile = sprintf('%s/sm0-swf-omni.flac',flacdir);

sess = '/homes/4/greve/links/sg1/xval/dng';
flacdir = '/homes/4/greve/links/sg1/xval/flac';
%flacfile = sprintf('%s/samc.flac',flacdir);


flaclist = '';

if(0)
flaclist = 'samc';
flaclist = strvcat(flaclist,'samc-rsyn');
flaclist = strvcat(flaclist,'samc-rsyn-f025');
flaclist = strvcat(flaclist,'samc-rsyn-f050');
flaclist = strvcat(flaclist,'samc-rsyn-f075');
flaclist = strvcat(flaclist,'samc-rsyn-f100');
end

flaclist = strvcat(flaclist,'samc-rsyn-snc-00-gn');
flaclist = strvcat(flaclist,'samc-rsyn-snc-05-gn');
flaclist = strvcat(flaclist,'samc-rsyn-snc-10-gn');

if(0)
flaclist = strvcat(flaclist,'samc-rsyn-snc-15');
flaclist = strvcat(flaclist,'samc-rsyn-snc-20');
flaclist = strvcat(flaclist,'samc-rsyn-snc-25');
flaclist = strvcat(flaclist,'samc-rsyn-snc-30');
flaclist = strvcat(flaclist,'samc-rsyn-snc-60');
end

nflacs = size(flaclist,1);

for nthflac = 1:nflacs
  flacname = deblank(flaclist(nthflac,:)); 
  fprintf('flac %s\n',flacname);
  
  %flacfile = sprintf('%s/samc-rsyn-sm5.flac',flacdir);
  flacfile = sprintf('%s/%s.flac',flacdir,flacname);
    
    
flac = fast_ldflac(flacfile);
flac.sess = sess;

flac.nthrun = 1;
flac = flac_customize(flac);
flac = flac_desmat(flac);

nruns = size(flac.runlist,1);
ncon = length(flac.con);

flarfxdir = sprintf('%s/%s/fla/%s/rfx',flac.sess,flac.fsd, ...
		    flac.name);
mkdirpcmd = sprintf('mkdir -p %s',flarfxdir);
unix(mkdirpcmd);

betamnsaved = 0;

for nthcon = 1:ncon

  C = flac.con(nthcon).C;
  J = size(C,1);
  fprintf('%d %s, J=%d\n',nthcon,flac.con(nthcon).name,J);
  if(J ~= 1) continue; end

  
  betasum = 0;
  gam = [];
  for nthrun = 1:nruns
    fprintf('nthrun = %d/%d\n',nthrun,nruns);
    flac.nthrun = nthrun;
    flac = flac_customize(flac);
    betarun = MRIread(flac.betafspec);
    if(~betamnsaved) betasum = betasum + betarun.vol; end
    betarun.vol = fast_vol2mat(betarun.vol);
    gamrun = C*betarun.vol;
    gam = [gam; gamrun];
  end % run

  X = ones(nruns,1);
  [beta rvar vdof] = fast_glmfitw(gam,X);
  [F dof1 dof2 ces] = fast_fratiow(beta,X,rvar,1);
  p = FTest(dof1, dof2, F);
  sig = -log10(p);
  sigmri = betarun;
  sigmri.vol = fast_mat2vol(sig,sigmri.volsize);

  condir = sprintf('%s/%s',flarfxdir,flac.con(nthcon).name);
  mkdirpcmd = sprintf('mkdir -p %s',condir);
  unix(mkdirpcmd);
  outfspec = sprintf('%s/fsig.mgh',condir);
  MRIwrite(sigmri,outfspec);
  
  Fmri = betarun;
  Fmri.vol = fast_mat2vol(F,Fmri.volsize); 
  outfspec = sprintf('%s/f.mgh',condir);
  MRIwrite(Fmri,outfspec);

  if(~betamnsaved)
    betamn = betarun;
    betamn.vol = betasum/nruns;
    betafspec = sprintf('%s/beta.mgh',flarfxdir);
    MRIwrite(betamn,betafspec);
    betamnsaved = 1;
  end
  
end % contrast

end % flac

fprintf('flacrfx done\n');



