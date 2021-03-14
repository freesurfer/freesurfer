function r = fast_bfileconvert(varargin)
% r = fast_bfileconvert(varargin)
% Converts a bfile into another bfile. Eg, bfloat into bshort,
% Little endian into big endian, etc.


%
% fast_bfileconvert.m
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

version = 'fast_bfileconvert.m @FS_VERSION@';
fprintf(1,'%s\n',version);
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

if(s.verbose)
  bfc_print_struct(s);
end

% ---------- Go through each slice -------------- %
nthslice = 1;
for slice = s.firstslice:s.lastslice
  fprintf(1,'%2d ',slice);
  if(rem(slice+1,10)==0) fprintf(1,'\n'); end

  %% Load the data %%
  %fname = sprintf('%s_%03d.%s',s.invol,slice,s.involext);
  %yi = fmri_ldbfile(fname);
  [yi mristruct] = fast_ldbslice(s.invol,slice);
  yi = yi([s.firstrow:s.lastrow]+1,...
          [s.firstcol:s.lastcol]+1,...
          [s.firstplane:s.lastplane]+1);
  if(~isempty(mristruct))
    mristruct.voldim(3) = s.lastslice-s.firstslice+1;
  end
  
  if(s.ln2log10) 
    p = sign(yi) .* exp(-abs(yi));
    iz = find(p==0); % do not take log of zero
    p(iz) = 1;
    yi = -sign(p).*log10(abs(p));
    yi(iz) = max(reshape1d(yi));
  end

  if(s.oddplanes) yi = yi(:,:,1:2:end);  end

  %% Save the data %%
  %outfname = sprintf('%s_%03d.%s',s.outvol,nthslice-1,s.outvolext);
  %fmri_svbfile(yi,outfname,s.outvolendian);
  fast_svbslice(yi,s.outvol,nthslice-1,s.outvolext,mristruct);

  nthslice = nthslice + 1;
end

if(s.rescaleshort)
  fprintf('Rescaling to fit in a bshort\n');
  y = fast_ldbslice(s.outvol);
  ymin = min(reshape1d(y));
  ymax = max(reshape1d(y));
  fprintf('ymin = %g, ymax = %g\n',ymin,ymax);  
  f = (2^15-1)/(ymax-ymin);
  y = f*(y-ymin);
  fast_svbslice(y,s.outvol,-1,s.outvolext,mristruct);  
end

r = 0;
fprintf(1,'fast_bfileconvert: completed SUCCESSFULLY\n');

return;
%----------------------------------------------------------%
%----------------------------------------------------------%
%----------------------------------------------------------%
%----------------------------------------------------------%


%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = bfc_struct;
  inputargs = varargin{1};
  ninputargs = length(inputargs);

  narg = 1;
  while(narg <= ninputargs)

    flag = deblank(inputargs{narg});
    narg = narg + 1;
    %fprintf(1,'Argument: %s\n',flag);
    if(~isstr(flag))
      flag
      fprintf(1,'ERROR: All Arguments must be a string\n');
      error;
    end

    switch(flag)

      case '-i',
        arg1check(flag,narg,ninputargs);
        s.invol = inputargs{narg};
        narg = narg + 1;

      case {'-firstrow', '-fr'}
        arg1check(flag,narg,ninputargs);
        s.firstrow = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-nrows', '-nr'}
        arg1check(flag,narg,ninputargs);
        s.nrows = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-firstcol', '-fc'}
        arg1check(flag,narg,ninputargs);
        s.firstcol = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-ncols', '-nc'}
        arg1check(flag,narg,ninputargs);
        s.ncols = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-firstslice', '-fs'}
        arg1check(flag,narg,ninputargs);
        s.firstslice = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-nslices', '-ns'}
        arg1check(flag,narg,ninputargs);
        s.nslices = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-firstplane', '-fp'}
        arg1check(flag,narg,ninputargs);
        s.firstplane = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-nplanes', '-np'}
        arg1check(flag,narg,ninputargs);
        s.nplanes = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-oddframes','-oddplanes'}
        s.oddplanes = 1;

      case '-o',
        arg1check(flag,narg,ninputargs);
        s.outvol = inputargs{narg};
        narg = narg + 1;

      case '-oext',
        arg1check(flag,narg,ninputargs);
        s.outvolext = inputargs{narg};
        narg = narg + 1;

      case '-oend',
        arg1check(flag,narg,ninputargs);
        s.outvolendian = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case '-ln2log10', 
        s.ln2log10 = 1;

      case '-rescale-short', 
        s.rescaleshort = 1; % rescale to fit in bshort

      case '-monly', % ignore
        arg1check(flag,narg,ninputargs);
        narg = narg + 1;

      case '-verbose',
        s.verbose = 1;

      case {'-debug','-echo','-umask'}, % ignore

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
function arg1check(flag,nflag,nmax)
  if(nflag>nmax) 
    fprintf(1,'ERROR: Flag %s needs one argument',flag);
    error;
  end
return;

%--------------------------------------------------%
%% Print Usage 
function print_usage(dummy)

  fprintf(1,'USAGE:\n');
  fprintf(1,'  fast_bfileconvert\n');
  fprintf(1,'     -i invol \n');
  fprintf(1,'     -firstslice sliceno : 0 \n');
  fprintf(1,'     -nslices    nslices : auto \n');
  fprintf(1,'     -firstplane planeno : 0 \n');
  fprintf(1,'     -nplanes    nplanes : auto \n');
  fprintf(1,'     -o    outvol \n');
  fprintf(1,'     -oext outvolext \n');
  fprintf(1,'     -oend outvolendianness \n');

return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = bfc_struct
  s.invol     = '';
  s.involext     = '';
  s.involendian  = -1;
  s.firstslice = 0;
  s.nrows    = -1;
  s.firstrow = 0;
  s.ncols    = -1;
  s.firstcol = 0;
  s.nslices    = -1;
  s.firstplane = 0;
  s.nplanes    = -1;
  s.oddplanes  = 0;
  s.outvol     = '';
  s.outvolext     = '';
  s.outvolendian  = -1;
  s.verbose = 0;
  s.ln2log10 = 0;
  s.rescaleshort = 0;
return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Print data structure
function s = bfc_print_struct(s,fid)
  if(nargin == 1) fid = 1; end

  fprintf(fid,'invol        %s\n',s.invol);
  fprintf(fid,'involext     %s\n',s.involext);
  fprintf(fid,'involendian  %d\n',s.involendian);
  fprintf(fid,'firstslice   %d\n',s.firstslice);
  fprintf(fid,'nslices      %d\n',s.nslices);
  fprintf(fid,'firstplane   %d\n',s.firstplane);
  fprintf(fid,'nplanes      %d\n',s.nplanes);
  fprintf(fid,'oddplanes    %d\n',s.oddplanes);
  fprintf(fid,'outvol       %s\n',s.invol);
  fprintf(fid,'outvolext    %s\n',s.involext);
  fprintf(fid,'outvolendian %d\n',s.involendian);

return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Check argument consistency, etc %%%
function s = check_params(s)

  fprintf(1,'Checking Parameters\n');
  if(size(s.invol,1) < 1) 
    msg = sprintf('ERROR: must have an input volume');
    qoe(msg);error(msg);
  end

  [nrows ncols ntp fs ns s.involendian s.involext] = fmri_bfiledim(s.invol);
  if(s.nrows == -1)        s.nrows = nrows; end
  if(s.ncols == -1)        s.ncols = ncols; end
  if(s.nslices == -1)      s.nslices = ns; end
  if(s.nplanes == -1)      s.nplanes = ntp; end
  if(isempty(s.outvol))    s.outvol       = s.invol; end
  if(isempty(s.outvolext)) s.outvolext    = s.involext;  end
  if(s.outvolendian == -1) s.outvolendian = s.involendian; end

  s.lastrow   = s.firstrow + s.nrows - 1;
  s.lastcol   = s.firstcol + s.ncols - 1;
  s.lastslice = s.firstslice + s.nslices - 1;
  s.lastplane = s.firstplane + s.nplanes - 1;

  if(s.nslices < 1) 
     msg = sprintf('ERROR: nslices = %d, must be > 0',s.nslices);
     qoe(msg);error(msg);
  end
  if(s.firstslice < 0) 
     msg = sprintf('ERROR: firstslice (%d) < 0',s.firstslice);
     qoe(msg);error(msg);
  end
  if(s.lastslice >= ns) 
     msg = sprintf('ERROR: last slices (%d) >= nslices avail (%d)',...
                   s.lastslice,ns);
     qoe(msg);error(msg);
  end

  if(s.nplanes < 1) 
     msg = sprintf('ERROR: nplanes = %d, must be > 0',s.nplanes);
     qoe(msg);error(msg);
  end
  if(s.firstplane < 0) 
     msg = sprintf('ERROR: firstplane (%d) < 0',s.firstplane);
     qoe(msg);error(msg);
  end
  if(s.lastplane >= ntp) 
     msg = sprintf('ERROR: last plane (%d) >= nplanes avail (%d)',...
                   s.lastplane,ntp);
     qoe(msg);error(msg);
  end

%  if(strcmp(s.invol,s.outvol) & strcmp(s.involext,s.outvolext) & ...
%     s.involendian == s.outvolendian )
%     msg = sprintf('INFO: nothing to do.');
%     qoe(msg);error(msg);
%  end
     
  if(strcmp(s.invol,s.outvol))
    if(s.verbose) fprintf(1,'INFO: overwriting input volume\n'); end
  end

return;


