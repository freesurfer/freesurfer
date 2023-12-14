% fcseedcor.m

if(0)
sesslist = char(textread('sidlist','%s'));
fsd = 'rest';
outfile = 'seedcor.dat';
seedlist = strvcat('fc.licc.dat','fc.postc.rh.dat','fc.calc.rh.dat','fc.cblum.dat');
xreglist = strvcat('global.waveform.dat','fc.vcsf.dat','fc.wm.dat','mcprextreg');
nperxreg = [1 1 1 6];
HPFCutOffHz = .01;
LPFCutOffHz = .20;
nskip = 4;
%------------------------------------------------------------------------------------
LPFOrderRatio = 10;
nsim = 0;
%------------------------------------------------------------------------------------
end



nsess = size(sesslist,1);
for nthsess = 1:nsess
  sess = deblank(sesslist(nthsess,:));
  fprintf('%2d %s\n',nthsess,sess);

fsdpath = sprintf('%s/%s',sess,fsd);
runlist = fast_runlist(fsdpath);
nruns = size(runlist,1);
nseeds = size(seedlist,1);
nxregs = [];

cmtxRuns = zeros(nseeds,nseeds,nruns);
pcmtxRuns = zeros(nseeds,nseeds,nruns);
for nthrun = 1:nruns

  xregs = [];
  for nthxreg = 1:size(xreglist,1)
    xregname = deblank(xreglist(nthxreg,:));
    fname = sprintf('%s/%s/%s',fsdpath,runlist(nthrun,:),xregname);
    s = load(fname,'-ascii');
    if(isempty(s)) return; end
    xreg = s(nskip+1:end,1:nperxreg(nthxreg));
    xregmn = mean(xreg);
    xreg = xreg - repmat(xregmn,[size(xreg,1) 1]);
    xregs = [xregs xreg];
  end

  seeds = [];
  for nthseed = 1:nseeds
    seedname = deblank(seedlist(nthseed,:));
    fname = sprintf('%s/%s/%s',fsdpath,runlist(nthrun,:),seedname);
    s = load(fname,'-ascii');
    if(isempty(s)) return; end
    seeds(:,nthseed) = s(nskip+1:end,1); % Only take the 1st
  end

  ntp = size(seeds,1);

  if(nthrun == 1)
    % Get TR so that HPF Cutoff can be computed
    fname = sprintf('%s/%s/f',fsdpath,runlist(nthrun,:));
    f = MRIread(fname,1);
    if(isempty(f)) return; end
    TR = f.tr/1000; % TR in sec
  end
  polyorder = fast_polyorder(ntp,TR,HPFCutOffHz);
  %fprintf('Run %d, polyorder = %d\n',nthrun,polyorder);
  Xp = fast_polytrendmtx(1,ntp,1,polyorder);
  Xr = [xregs Xp];

  LPFOrder = round(ntp/LPFOrderRatio);
  if(LPFCutOffHz > 0) 
    LPF = fast_lpfmtx(LPFCutOffHz,TR,ntp,LPFOrder);
  end

  % Compute partial correlation coefficient
  for nthseed = 1:nseeds-1
    for mthseed = nthseed+1:nseeds
      y = seeds(:,nthseed);
      X = [seeds(:,mthseed) Xr];
      if(LPFCutOffHz > 0) 
	y = LPF*y;
	X = LPF*X;
      end
      [beta rvar] = fast_glmfit(y,X);
      C = zeros(1,size(X,2));
      C(1) = 1;
      rho = fast_glm_pcc(beta,X,C,rvar);
      cmtxRuns(nthseed,mthseed,nthrun) = rho;
      cmtxRuns(mthseed,nthseed,nthrun) = rho;
      
      if(nsim > 0)
	rhosim = zeros(nsim,1);
	for nthsim = 1:nsim
	  y = randn(size(seeds(:,nthseed)));
	  if(LPFCutOffHz > 0) 
	    y = LPF*y;
	  end
	  [beta rvar] = fast_glmfit(y,X);
	  C = zeros(1,size(X,2));
	  C(1) = 1;
	  rhosim(nthsim) = fast_glm_pcc(beta,X,C,rvar);
	end
	pcmtxRuns(nthseed,mthseed,nthrun) = length(find(rhosim > abs(rho)))/nsim;
	pcmtxRuns(mthseed,nthseed,nthrun) = length(find(rhosim > abs(rho)))/nsim;
      end      
    end
  end

  fname = sprintf('%s/%s/%s.log',fsdpath,runlist(nthrun,:),outfile);
  fp = fopen(fname,'w');
  fprintf(fp,'Log file for fcseedcor.m\n');
  fprintf(fp,'Date %s\n',date);
  fprintf(fp,'TR %g\n',TR);
  fprintf(fp,'ntp %g\n',ntp);
  fprintf(fp,'HPFCutOffHz %g\n',HPFCutOffHz);
  fprintf(fp,'PolyOrder %d\n',polyorder);
  fprintf(fp,'LPFCutOffHz %g\n',LPFCutOffHz);
  fprintf(fp,'LPFOrderRatio %d\n',LPFOrderRatio);
  fprintf(fp,'LPFOrder %d\n',LPFOrder);
  fprintf(fp,'nskip %d\n',nskip);
  fprintf(fp,'nsim %d\n',nsim);
  fprintf(fp,'\n');
  for nthxreg = 1:size(xreglist,1)
    xregname = deblank(xreglist(nthxreg,:));
    fprintf(fp,'%2d xreg %s %2d\n',nthxreg,xregname,nperxreg(nthxreg));
  end
  for nthseed = 1:size(seedlist,1)
    seedname = deblank(seedlist(nthseed,:));
    fprintf(fp,'%2d seed %s\n',nthseed,seedname);
  end
  fclose(fp);

  fname = sprintf('%s/%s/%s',fsdpath,runlist(nthrun,:),outfile);
  fp = fopen(fname,'w');
  for nthseed = 1:nseeds
    for mthseed = 1:nseeds
      fprintf(fp,'%8.5f ',cmtxRuns(nthseed,mthseed,nthrun));
    end
    fprintf(fp,'\n');
  end
  fclose(fp);

  if(nsim > 0) 
    fname = sprintf('%s/%s/p%s',fsdpath,runlist(nthrun,:),outfile);
    fp = fopen(fname,'w');
    for nthseed = 1:nseeds
      for mthseed = 1:nseeds
	if(nthseed == mthseed) 	fprintf(fp,'%7.5f ',1);
	else fprintf(fp,'%7.5f ',pcmtxRuns(nthseed,mthseed,nthrun));
	end
      end
      fprintf(fp,'\n');
    end
    fclose(fp);
  end
  
end

cmtx = mean(cmtxRuns,3);

end

fprintf('fast_fcseedcor done\n');
if(~monly) quit; end
return; 







