function r = fast_intergroupavg(varargin)
% Name: fast_intergroupavg
% Purpose: implements inter-group averaging
%
% Author: Douglas Greve
% Questions or Comments: analysis-bugs@nmr.mgh.harvard.edu
% r = fast_isxavg_fe(varargin)


%
% fast_intergroupavg.m
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

r = 1;

%% Print useage if there are no arguments %%
if(nargin == 0)
  print_usage;
  return;
end

%% Parse the arguments %%
s = parse_args(varargin);
if(isempty(s)) return; end
s = check_params(s);
if(isempty(s)) return; end

if(s.synth < 0) s.synth = sum(100*clock); end
fprintf('SynthSeed = %10d\n',s.synth);

fprintf(1,'_______________________________________________\n');
fprintf(1,'Intergroup Averaging Parameters\n');
iga_print_struct(s,1);
fprintf(1,'_______________________________________________\n');

[ns nr nc nf] = fmri_bvoldim(s.avg1id);
nv = nr*nc;

lastslice = s.firstslice + s.nslices - 1;

for slice = s.firstslice:lastslice

  % fprintf('Processing Slice %d\n',slice);
  fprintf('%2d ',slice);

  fname = sprintf('%s_%03d.bfloat',s.avg1id,slice);
  avg1 = fmri_ldbfile(fname);
  avg1 = reshape(avg1, [nv nf])';%'

  fname = sprintf('%s_%03d.bfloat',s.std1id,slice);
  std1 = fmri_ldbfile(fname);
  std1 = reshape(std1, [nv nf])';%'

  fname = sprintf('%s_%03d.bfloat',s.avg2id,slice);
  avg2 = fmri_ldbfile(fname);
  avg2 = reshape(avg2, [nv nf])';%'

  fname = sprintf('%s_%03d.bfloat',s.std2id,slice);
  std2 = fmri_ldbfile(fname);
  std2 = reshape(std2, [nv nf])';%'

  if(s.polarity > 0)
    ind = find(avg1 < 0);
    avg1(ind) = 1e-10;
    std1(ind) = 1e+10;
    ind = find(avg2 < 0);
    avg2(ind) = 1e-10;
    std2(ind) = 1e+10;
  end
  if(s.polarity < 0)
    ind = find(avg1 > 0);
    avg1(ind) = 1e-10;
    std1(ind) = 1e+10;
    ind = find(avg2 > 0);
    avg2(ind) = 1e-10;
    std2(ind) = 1e+10;
  end

  % Polarity2 - only set the average to 0; leave the std as
  % it was.
  if(s.polarity2 > 0)
    ind = find(avg1 < 0);
    fprintf('INFO: found %d voxels less than zero in group1\n',...
	    length(ind));
    avg1(ind) = 0;
    ind = find(avg2 < 0);
    fprintf('INFO: found %d voxels less than zero in group2\n',...
	    length(ind));
    avg2(ind) = 0;
  end
  if(s.polarity2 < 0)
    ind = find(avg1 > 0);
    fprintf('INFO: found %d voxels greater than zero in group1\n',...
	    length(ind));
    avg1(ind) = 0;
    ind = find(avg2 > 0);
    fprintf('INFO: found %d voxels greater than zero in group2\n',...
	    length(ind));
    avg2(ind) = 0;
  end

  stderr = sqrt( (std1.^2)/s.dof1 + (std2.^2)/s.dof2 );
  avgdiff = avg1-avg2;

  % Find those places where the stddev and the average are
  % zero. This can happen when the FOVs of the individuals
  % dont completely overlap with those of the registration 
  % subject. The avgdiff can equal zero if polarity is used
  % because two voxels could be less than zero in both
  % avg1 and avg2 in which case they would be set to the
  % same value. So the number of zeros in avgdiff may be 
  % much more than that in stderr.
  iz = find(avgdiff == 0);
  avgdiff(iz) = .0000001;
  iz = find(stderr == 0);
  stderr(iz) = 100000000;

  t = avgdiff ./ stderr;
  dof = s.dof1 + s.dof2 - 2;
  tsig = tTest(dof,reshape1d(t));
  tsig = reshape(tsig, size(t));
  tsig = sign(t).*tsig;

  % Find those voxels whose t value is so large that the significance
  % is quantized to  zero. Replace zero with a very small number. This
  % keeps the log10 from failing later on.
  iz = find(tsig==0); 
  tsig(iz) = 1e-100;

  if(~isempty(s.tid))
    fname = sprintf('%s_%03d.bfloat',s.tid,slice);
    tmp = reshape(t', [nr nc nf]); %'
    fmri_svbfile(tmp,fname);
  end

  if(~isempty(s.tsigid))
    fname = sprintf('%s_%03d.bfloat',s.tsigid,slice);
    tmp = -sign(tsig) .* log10(-abs(tsig));
    tmp = reshape(tmp', [nr nc nf]); %'
    fmri_svbfile(tmp,fname);
  end

  if(~isempty(s.tminsigid))
    fname = sprintf('%s_%03d.bfloat',s.tminsigid,slice);
    [tminsig itminsig] = min(abs(tsig),[],1);
    ind = sub2ind(size(tsig),itminsig,1:nv);
    tminsig = tsig(ind); % signed pmin
    tminsig = tminsig * nf; % Bonferoni Correction
    tmp = -sign(tminsig) .* log10(-abs(tminsig));
    tmp = reshape(tmp', [nr nc]); %'
    fmri_svbfile(tmp,fname);
  end

  clear tmp;

end % Loop over slices %

fprintf('\n');

fprintf(1,'fast_intergroupavg Completed \n');

return;

%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = iga_struct;;
  inputargs = varargin{1};
  ninputargs = length(inputargs);

  narg = 1;
  while(narg <= ninputargs)

    flag = deblank(inputargs{narg});
    narg = narg + 1;
    % fprintf(1,'Argument: %s\n',flag);
    if(~isstr(flag))
      flag
      fprintf(1,'ERROR: All Arguments must be a string\n');
      s = []; return;
    end

    switch(flag)

      case {'-firstslice', '-fs'}
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.firstslice = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-nslices', '-ns'}
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.nslices = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-avg1'},
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.avg1id = inputargs{narg};
        narg = narg + 1;

      case {'-std1'},
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.std1id = inputargs{narg};
        narg = narg + 1;

      case {'-dof1'}
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.dof1 = sscanf(inputargs{narg},'%d',1);
        if(s.dof1 < 1)
          fprintf(1,'ERROR: dof1 = %d, must be > 0\n',s.dof1);
          s = []; return;
        end
        narg = narg + 1;

      case {'-avg2'},
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.avg2id = inputargs{narg};
        narg = narg + 1;

      case {'-std2'},
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.std2id = inputargs{narg};
        narg = narg + 1;

      case {'-dof2'}
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.dof2 = sscanf(inputargs{narg},'%d',1);
        if(s.dof2 < 1)
          fprintf(1,'ERROR: dof2 = %d, must be > 0\n',s.dof2);
          s = []; return;
        end
        narg = narg + 1;

      case {'-polarity'}
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.polarity = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-polarity2'}
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.polarity2 = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-t'},
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.tid = inputargs{narg};
        narg = narg + 1;

      case {'-tsig'},
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.tsigid = inputargs{narg};
        narg = narg + 1;

      case {'-tminsig'},
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.tminsigid = inputargs{narg};
        narg = narg + 1;

      case '-synth',
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        s.synth  = sscanf(inputargs{narg},'%f',1);
        narg = narg + 1;

      case {'-monly', '-umask'}, % ignore
        if(arg1check(flag,narg,ninputargs)) s = []; return; end
        narg = narg + 1;

      otherwise
        fprintf(2,'ERROR: Flag %s unrecognized\n',flag);
        s = [];
        return;

    end % --- switch(flag) ----- %

  end % while(narg <= ninputargs)

return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Check that there is at least one more argument %%
function r = arg1check(flag,nflag,nmax)
  r = 0;
  if(nflag>nmax) 
    fprintf(1,'ERROR: Flag %s needs one argument',flag);
    r = 1;
  end
return;

%--------------------------------------------------%
%% Check that there are at least two more arguments %%
function r = arg2check(flag,nflag,nmax)
  r = 0;
  if(nflag+1 > nmax) 
    fprintf(1,'ERROR: Flag %s needs two arguments',flag);
    r = 1;
  end
return;

%--------------------------------------------------%
%% Check argument consistency, etc %%%
function s = check_params(s)

  fprintf(1,'Checking Parameters\n');

  if(isempty(s.avg1id))
    fprintf(1,'ERROR: must specify avg1\n');
    s = []; return;
  end

  if(isempty(s.std1id))
    fprintf(1,'ERROR: must specify std1\n');
    s = []; return;
  end

  if(isempty(s.avg2id))
    fprintf(1,'ERROR: must specify avg2\n');
    s = []; return;
  end

  if(isempty(s.std2id))
    fprintf(1,'ERROR: must specify std2\n');
    s = []; return;
  end

  if(isempty(s.tid) & isempty(s.tsigid) & isempty(s.tminsigid))
    fprintf(1,'ERROR: must specify an output (t, tsig, tminsig)\n');
    s = []; return;
  end

  if(s.dof1 == -1)
    fname = sprintf('%s.dof',s.std1id);
    fid = fopen(fname,'r');
    if(fid == -1)
      fprintf(1,'ERROR: cannot open %s\n',fname);
      s = []; return;
    end
    s.dof1 = fscanf(fid,'%d',1);
    fclose(fid);      
  end
  if(s.dof1 < 1)
    fprintf(1,'ERROR: dof1 = %d, must be >= 0\n',s.dof1);
    s = []; return;
  end

  if(s.dof2 == -1)
    fname = sprintf('%s.dof',s.std2id);
    fid = fopen(fname,'r');
    if(fid == -1)
      fprintf(1,'ERROR: cannot open %s\n',fname);
      s = []; return;
    end
    s.dof2 = fscanf(fid,'%d',1);
    fclose(fid);      
  end
  if(s.dof2 < 1)
    fprintf(1,'ERROR: dof2 = %d, must be >= 0\n',s.dof2);
    s = []; return;
  end

  if(s.dof1 + s.dof2 < 3)
    fprintf(1,'ERROR: dof1+dof2 = %d, must be >= 3\n',s.dof1+s.dof2);
    s = []; return;
  end
    


  [a1ns a1nr a1nc a1nf] = fmri_bvoldim(s.avg1id);
  [s1ns s1nr s1nc s1nf] = fmri_bvoldim(s.std1id);
  [a2ns a2nr a2nc a2nf] = fmri_bvoldim(s.avg2id);
  [s2ns s2nr s2nc s2nf] = fmri_bvoldim(s.std2id);

  if(a1ns ~= a2ns)
    fprintf(1,'ERROR: avg1 nslices != avg2 nslices (%d,%d)\n',a1ns,a2ns);
    s = []; return;
  end

  if(a1nf ~= a2nf)
    fprintf(1,'ERROR: avg1 nframes != avg2 frames (%d,%d)\n',a1nf,a2nf);
    s = []; return;
  end

  if(s.nslices < 1) s.nslices = a1ns; end


return;

%--------------------------------------------------%
%% Print Usage 
function print_usage
  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_intergroupavg \n');
  fprintf(1,'     -avg1  avg1id \n');
  fprintf(1,'     -std1  std1id \n');
  fprintf(1,'     -dof1  dof1 \n');
  fprintf(1,'     -avg2  avg2id \n');
  fprintf(1,'     -std2  std2id \n');
  fprintf(1,'     -dof2  dof2 \n');

  fprintf(1,'     -t       tid\n');
  fprintf(1,'     -tsig    tsigid\n');
  fprintf(1,'     -tminsig tminsigid\n');

  fprintf(1,'     -firstslice sliceno  \n');
  fprintf(1,'     -nslices    nslices  \n');
return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = iga_struct;
  s.avg1id    = '';
  s.std1id    = '';
  s.dof1      = -1;
  s.avg2id    = '';
  s.std2id    = '';
  s.dof2      = -1;
  s.firstslice = 0;
  s.nslices    = -1;
  s.tid = '';
  s.tsigid = '';
  s.tminsigid = '';
  s.polarity  = 0;
  s.synth   = 0;
return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Print data structure
function s = iga_print_struct(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'avg1id      %s\n',s.avg1id);
  fprintf(fid,'std1id      %s\n',s.std1id);
  fprintf(fid,'dof1        %d\n',s.dof1);

  fprintf(fid,'avg2id      %s\n',s.avg2id);
  fprintf(fid,'std2id      %s\n',s.std2id);
  fprintf(fid,'dof2        %d\n',s.dof2);

  fprintf(fid,'tid         %s\n',s.tid);
  fprintf(fid,'tsigid      %s\n',s.tsigid);
  fprintf(fid,'tminsigid   %s\n',s.tminsigid);

  fprintf(fid,'synth       %s\n',s.synth);

return;
%--------------------------------------------------%



