% flacffx.m
% $Id: flacffx.m,v 1.1 2005/01/11 19:00:24 greve Exp $
tic;

%flacfile = '/homes/4/greve/fbirn-hp-fsfast/flac/bh.flac';
%sess = '/homes/4/greve/fbirn-hp-fsfast/mgh-data/mgh-103.1';
%sess = '/homes/4/greve/fbirn-hp-fsfast/mgh-data/mgh-101.1';
%flacfile = '/homes/4/greve/links/amn/flac/fn.flac';
%sess = '/homes/4/greve/links/amn/AMN_12';
%flacfile = '/homes/4/greve/fbirn-hp-fsfast/flac/sm-wgn.flac';

%sess = '/homes/4/greve/links/sg1/skb/skb-101.1';
%flacdir = '/homes/4/greve/links/sg1/skb/flac/';
%flacfile = sprintf('%s/fn-swf-nvr.flac',flacdir);

%sess = '/homes/4/greve/links/sg1/dng081403/dng-siemens';
%flacdir = '/homes/4/greve/links/sg1/dng081403/flac';
%flacfile = sprintf('%s/sm0-swf-omni.flac',flacdir);

sess = '/homes/4/greve/links/sg1/xval/dng';
flacdir = '/homes/4/greve/links/sg1/xval/flac';
%flacfile = sprintf('%s/samc.flac',flacdir);
flacfile = sprintf('%s/sem_assoc-rsyn-swfnbhd.flac',flacdir);
    
    
flac = fast_ldflac(flacfile);
flac.sess = sess;

flac.nthrun = 1;
flac = flac_customize(flac);
flac = flac_desmat(flac);

nruns = size(flac.runlist,1);
ncon = length(flac.con);

flaffxdir = sprintf('%s/%s/fla/%s/ffx',flac.sess,flac.fsd, ...
		    flac.name);
mkdirpcmd = sprintf('mkdir -p %s',flaffxdir);
unix(mkdirpcmd);

dof = 0;
for nthrun = 1:nruns
  fprintf('nthrun = %d/%d\n',nthrun,nruns);
  flac.nthrun = nthrun;
  flac = flac_customize(flac);
  flac = flac_desmat(flac);
  dofrun = size(flac.X,1) - size(flac.X,2);
  dof = dof + dofrun;
  
  fladir = sprintf('%s/%s/fla/%s/%s',...
		   flac.sess,flac.fsd,flac.name,flac.runlist(nthrun,:));
  
  betafspec = sprintf('%s/beta',fladir);
  betarun = MRIread(betafspec);
  if(nthrun == 1)
    betamn = betarun;
    betamn.vol = 0;
  end
  betamn.vol  = betamn.vol + betarun.vol/nruns;
  betarun.vol = fast_vol2mat(betarun.vol);
  
  rvarfspec = sprintf('%s/rvar',fladir);
  rvarrun = MRIread(rvarfspec);
  ssrrun = dofrun*rvarrun.vol;
  if(nthrun == 1)
    ssr = rvarrun;
    ssr.vol = 0;
  end
  ssr.vol = ssr.vol + ssrrun;

  matfile = sprintf('%s/flac.mat',fladir);
  flacproc = load(matfile);
  nseg = size(flacproc.nacfseg,2);
  
  for nthcon = 1:ncon
    C = flac.con(nthcon).C;
    J = size(C,1);
    gamrun = C*betarun.vol;
    if(nthrun == 1)  
      gam(nthcon) = betarun; 
      gam(nthcon).vol = 0;
    end
    gam(nthcon).vol = gam(nthcon).vol + gamrun;
    
    for nthseg = 1:nseg+1
      if(nthseg == 1) S = eye(flac.ntp);
      else
	nacf = flacproc.nacfseg(:,nthseg-1);
	S = toeplitz(nacf);
      end
      cixtxctrun = C*inv(flac.X'*inv(S)*flac.X)*C';
      if(nthrun == 1) 
	flac.con(nthcon).cixtxct(:,:,nthseg) = zeros(J); 
      end
      flac.con(nthcon).cixtxct(:,:,nthseg) = ...
	  flac.con(nthcon).cixtxct(:,:,nthseg) + cixtxctrun;
    end % seg
  
  end % con
  
end % run

rvar = ssr;
rvar.vol = ssr.vol/dof;

outfspec = sprintf('%s/beta.%s',flaffxdir,flac.format);
MRIwrite(betamn,outfspec);

outfspec = sprintf('%s/rvar.%s',flaffxdir,flac.format);
MRIwrite(rvar,outfspec);


for nthcon = 1:ncon
  fprintf('con %s %d/%d\n',flac.con(nthcon).name,nthcon,ncon);
  C = flac.con(nthcon).C;
  J = size(C,1);
  F = gam(nthcon);
  F.vol = zeros(size(gam(nthcon).vol(1,:)));
  for nthseg = 1:nseg+1
    indseg = find(flacproc.acfseg.vol == nthseg-1);
    cixtxct = flac.con(nthcon).cixtxct(:,:,nthseg);
    tmp = (gam(nthcon).vol(:,indseg)' * inv(cixtxct))';
    if(J > 1)
      tmp = sum(tmp .* gam(nthcon).vol(:,indseg));
    else
      tmp = tmp .* gam(nthcon).vol(:,indseg);
    end      
    F.vol(indseg) = tmp./(J*rvar.vol(indseg)');
  end % seg
  p = F;
  p.vol = FTest(J, dof, F.vol);
  ind = find(p.vol == 0);
  p.vol(ind) = eps;
  sig = p;
  sig.vol = -log10(p.vol);

  condir = sprintf('%s/%s',flaffxdir,flac.con(nthcon).name);
  mkdirpcmd = sprintf('mkdir -p %s',condir);
  unix(mkdirpcmd);

  F.vol = fast_mat2vol(F.vol,F.volsize);
  outfspec = sprintf('%s/f.%s',condir,flac.format);
  MRIwrite(F,outfspec);
  
  sig.vol = fast_mat2vol(sig.vol,sig.volsize);
  outfspec = sprintf('%s/fsig.%s',condir,flac.format);
  MRIwrite(sig,outfspec);

end % con

fprintf('done\n');



