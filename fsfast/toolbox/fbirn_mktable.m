% fbirn_mktable.m - make data table for fbirn


%
% fbirn_mktable.m
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
pthresh = 3;

for partype = [3 2 1]

if(partype == 1) 
  partypeid = 'sm'; 
  runnolist = [3 5 7 10];
end
if(partype == 2) 
  partypeid = 'bh'; 
  runnolist = [6 8];
end
if(partype == 3) 
  partypeid = 'rest'; 
  runnolist = [4 9];
end

fprintf('----------------------------------------------\n');
fprintf('----------------------------------------------\n');
fprintf('----------------------------------------------\n');
  fprintf('partype = %s (%g)\n',partypeid,toc);

ana   = sprintf('%s-sm5-per',partypeid);
mcana = sprintf('%s-mc-per',partypeid);

topdir = '/autofs/space/greve_002/users/greve/fbirn-hp-fsfast';
tablefile = sprintf('%s/%s.table',topdir,partypeid);

% site subj visit partype run   nbrain %over %mag cnr %rstd ar1 pmc mcmax
fp = fopen(tablefile,'w');
fprintf(fp,['site subj visit partype run nbrain pctover magpct cnr' ...
	    ' rstdpct ar1 pmc mcmax\n']);
fclose(fp);

sitelist = '';
sitelist = strvcat(sitelist,'bwh');     % 1
sitelist = strvcat(sitelist,'dunc');    % 2
sitelist = strvcat(sitelist,'dunc4t');  % 3
sitelist = strvcat(sitelist,'iowa');    % 4
sitelist = strvcat(sitelist,'mgh');     % 5
sitelist = strvcat(sitelist,'min');     % 6
sitelist = strvcat(sitelist,'nm');      % 7
sitelist = strvcat(sitelist,'stan');    % 8
sitelist = strvcat(sitelist,'uci');     % 9
sitelist = strvcat(sitelist,'ucla');    % 10
sitelist = strvcat(sitelist,'ucsd');    % 11

nsites = size(sitelist,1);

for site = 1:nsites
  siteid = deblank(sitelist(site,:));

  if(site ~= 2 & site ~= 9) continue; end %%%%%%%%%%%%%%

  for subj = [1 3:6]

    if(subj ~= 3) continue; end %%%%%%%%%%%%%

    for visit = 1:2

      if(visit ~= 2) continue; end %%%%%%%%%%%%%

      if(strcmp(siteid,'mgh')==0)
	sessid = sprintf('%s-data/%s-10%d.%d',siteid,siteid,subj,visit);
      else
	sessid = sprintf('mgh-10%d.%d',subj,visit);
      end	

      fprintf('sessid %s  (%g)\n',sessid,toc);
      if(~direxists(sessid))
	fprintf(' ..... skipping\n');
	continue;
      end
      
      maskdir = sprintf('%s/bold/masks',sessid);
      if(~direxists(maskdir))
	fprintf('   no mask, skipping\n');
	nthrun = nthrun + 1;
	continue;
      end

      maskstem = sprintf('%s/brain',maskdir);
      mask = fast_ldbslice(maskstem);
      if(isempty(mask)) return; end
      mask(:,:,[1 end]) = 0; % exclude ends
      indmask_in  = find(mask==1);
      indmask_out = find(mask==0);
      nbrain = length(indmask_in);
      
      nthrun = 1;
      for runno = runnolist
	anadir = sprintf('%s/bold/%s-%03d',sessid,ana,runno);
	mcanadir = sprintf('%s/bold/%s-%03d',sessid,mcana,runno);

	if(~direxists(anadir) | ~direxists(mcanadir))
	  fprintf('   Skipping run %d\n',runno);
	  nthrun = nthrun + 1;
	  continue;
	end
	
	hoffsetstem = sprintf('%s/h-offset',anadir);
	hoffset = fast_ldbslice(hoffsetstem);
	if(isempty(hoffset)) return; end
	indz = find(hoffset==0);
	hoffset(indz) = 1e10;
	
	rvarstem = sprintf('%s/estsnr/rvar',anadir);
	rvar = fast_ldbslice(rvarstem);
	if(isempty(rvar)) return; end
	rvarpct = 100*rvar./hoffset;
	rstdpctmn = sqrt(mean(rvarpct(indmask_in)));

	ar1stem = sprintf('%s/estsnr/ar1',anadir);
	ar1 = fast_ldbslice(ar1stem);
	if(isempty(ar1)) return; end
	ar1mn = mean(ar1(indmask_in));

	pomnistem = sprintf('%s/omnibus/fsig',anadir);
	pomni = fast_ldbslice(pomnistem);
	if(isempty(pomni)) return; end
	pomni(indmask_out) = 0;
	indover = find(abs(pomni) > pthresh);
	nover = length(indover);
	pctover = 100*nover/nbrain;
	
	cnrstem = sprintf('%s/omnibus/f',anadir);
	cnr = fast_ldbslice(cnrstem);
	if(isempty(cnr)) return; end
	cnrmn = mean(cnr(indover));

	magpctstem = sprintf('%s/omnibus/magpct',anadir);
	magpct = fast_ldbslice(magpctstem);
	if(isempty(magpct)) return; end
	magpctmn = mean(magpct(indover));
	
	pmcomnistem = sprintf('%s/omnibus/fsig',mcanadir);
	pmcomni = fast_ldbslice(pmcomnistem);
	if(isempty(pmcomni)) return; end
	pmc = 10.^(-abs(pmcomni));
	pmcmin = min(pmc);
	pmcminbc = 1 - (1-pmcmin).^length(pmc); % bonferonni correction
	if(pmcminbc == 0) pmcminbc = pmcmin*length(pmc); end
	if(pmcminbc == 0) pmcminbc = eps; end
	pmcminlog10 = -log10(pmcminbc);

	mcdatfile = sprintf('%s/bold/%03d/fmc.mcdat',sessid,runno);
	mcdat = load(mcdatfile);
	if(isempty(mcdat))
	  fprintf('ERROR: loading %s\n',mcdatfile);
	  return;
	end
	mcmax = max(mcdat(:,10));

	% site subj visit partype run   nbrain %over %mag cnr %rstd ar1 pmc mcmax
	fmt = '%2d %d %d %d %2d  %5d %5.2f %5.2f %5.2f %4.2f %4.2f %4.1f %3.2f\n';
	dat = [site subj visit partype nthrun     nbrain pctover ...
	       magpctmn cnrmn rstdpctmn  ar1mn  pmcminlog10 mcmax];
	fp = fopen(tablefile,'a');
	fprintf(fp,fmt,dat);
	fclose(fp);
	fprintf(fmt,dat);

	% Make a local copy %
	fname = sprintf('%s/sum.table',anadir);
	fp = fopen(fname,'w');
	fprintf(fp,fmt,dat);
	fclose(fp);
	
	nthrun = nthrun + 1;
      end
      
    end
  end
end

end % partype


fprintf('fbirn_mktable: done %g\n',toc);
