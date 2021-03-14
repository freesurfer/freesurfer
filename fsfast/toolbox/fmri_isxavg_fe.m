% Name: fmri_isxavg_fe
% Purpose: implements fixed-effects intersession averaging
%          for output of selxavg
% Author: Douglas Greve
% Questions or Comments: analysis-bugs@nmr.mgh.harvard.edu


%
% fmri_isxavg_fe.m
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

%%%% These variables must be specified %%%%%%%%
% InStemList
% FirstSlice
% NSlices
% OutStem
% weight
% pctsigch
% truncsign

if( exist(deblank('InStemList')) ~= 1)
  msg = sprintf('Variable InStemList does not exist');
  qoe(msg);error(msg);
end

nsessions = size(InStemList,1);

LastSlice = FirstSlice + NSlices - 1;

for slice = FirstSlice:LastSlice

  fprintf('%3d ',slice);
  %fprintf('%2d ',slice);
  eVarSum = 0;
  DOFSum  = 0;
  hAvgSum = 0;
  hCovSum = 0;

  %fprintf('Session ');
  for session = 1:nsessions,
    %fprintf('%2d ',session);
    InStem  = deblank(InStemList(session,:));
    InSA    = sprintf('%s_%03d.bfloat',InStem,slice);
    DatFile = sprintf('%s.dat',InStem);
    InHOffset  = sprintf('%s-offset_%03d.bfloat',InStem,slice);
    
    [hAvg eVar hd] = fast_ldsxabfile(InSA);
    if(session == 1)
      hd0 = hd;
    else
      if(hd0.Nnnc ~= hd.Nnnc)
	fprintf('\n\nERROR: session %d has a different number of\n',session);
	fprintf('  conditions than session 1\n\n\n');
	return;
      end
      if(size(hAvg,1) ~= size(hAvgSum,1) | ...
	size(hAvg,2) ~= size(hAvgSum,2)  | ...
	size(hAvg,3) ~= size(hAvgSum,3) )
	fprintf(['ERROR: dimension mismatch between session '...
		 'number %d and previous sessions\n'],session);
	fprintf(['This usually happens when an analysis has been '...
		 'redefined but selxavg has not been re-run on all '...
		 'the subjects.\n']);
	return;
      end
    end
    %hAvg = randn(size(hAvg));
      
    % Truncate all values of the specified sign to 0
    if( ~isempty(truncsign) )
      if( strcmpi(truncsign,'pos') )
        ind = find(hAvg > 0);
      end
      if( strcmpi(truncsign,'neg') )
        ind = find(hAvg < 0);
      end
      %fprintf('ntrunc = %d\n',length(ind));
      if( ~isempty(ind) ) 
        hAvg(ind) = 0; 
	% What to do with eVar when truncating averages??
        %  Leave it - keep noise (OK for post-random FX because thown away)
        %  Truncate it - which ones, hAvg can have mult planes
      end
    end

    if(pctsigch)
      hoffset = fmri_ldbfile(InHOffset);
      ind = find(hoffset == 0);
      hoffset(ind) = 10^10;
      eVar = eVar ./(hoffset.^2);
      hofftmp = repmat(hoffset,[1 1 size(hAvg,3)]);
      hAvg = hAvg./hofftmp;
      clear hoffset hofftmp;
    end

    hCov = hd.hCovMtx;
    eVarSum = eVarSum + eVar * abs(weight(session)) * hd.DOF;
    clear eVar;
    hAvgSum = hAvgSum + hAvg*(weight(session) * hd.DOF);
    clear hAvg;
    hCovSum = hCovSum + inv(hCov) * abs(weight(session)) ;
    DOFSum  = DOFSum + hd.DOF * abs(weight(session));
  end % loop over sessions %
  % fprintf('  ---\n');

  DOFSum = round(DOFSum);

  hAvgGrp = hAvgSum/DOFSum;
  eVarGrp = eVarSum/DOFSum;
  hCovGrp = inv(hCovSum);

  hd.hCovMtx = hCovGrp;
  hd.DOF = DOFSum;
  [ySA dof] = fmri_sxa2sa(eVarGrp,hCovGrp,hAvgGrp,hd);

  OutSA = sprintf('%s_%03d.bfloat',OutStem,slice); 
  fmri_svbfile(ySA, OutSA); 

  OutDat = sprintf('%s.dat',OutStem);
  fmri_svdat2(OutDat,hd);

end % Loop over slices %
fprintf('\n');

fmri_touch(okfile);
fprintf(1,'fmri_isavg_fe completed \n\n');






