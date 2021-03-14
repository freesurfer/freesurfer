function r = fast_volview(varargin)
% Name: fast_volview(varargin)
% Author: Douglas Greve
% Questions or Comments: analysis-bugs@nmr.mgh.harvard.edu


%
% fast_volview.m
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

%% Print usage if there are no arguments %%
if(nargin == 0)  print_usage;  return; end

%% Parse the arguments %%
s = parse_args(varargin);
if(isempty(s)) return; end
%s = check_params(s);
%if(isempty(s)) return; end

if(s.redraw)
  fprintf('Redraw\n');
  s.displayimage = s.vol(:,:,s.curpos(3)+1,s.curpos(4)+1);
  s.himage = imagesc(s.displayimage);
  s.redraw = 0;
end

set(s.hfig,'userdata',s);

return;
%--------------------------------------------------%
%---- end main ------------------------------------%
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = main_struct
  s.hfig        = [];
  s.himage      = [];
  s.haxis       = [];
  s.vol         = [];
  s.volid       = [];
  s.volres      = [1 1 1 1];
  s.title       = 'untitled';
  s.mosaic      = 0;
  s.curpos      = [0 0 0 0];
  s.setpos      = [0 0 0 0];
  s.curchar     = '';
  s.displayimage = [];
  s.redraw       = 0;
  s.mousedown    = 0;
  s.motion       = 1;
return;
%--------------------------------------------------%

%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  s = main_struct;
  inputargs = varargin{1};
  ninputargs = length(inputargs);

  narg = 1;
  flag = deblank(inputargs{narg});
  if(strcmp(flag,'-h'))
    narg = narg + 1;
    if(argNcheck(flag,narg,ninputargs,1)) s = []; return; end
    hfig = inputargs{narg};
    if(isempty(hfig)) s = initfigure;
    else              s = get(hfig,'userdata');     
    end
    narg = narg + 1;
  else
    s = get(gcf,'userdata');     
  end

  while(narg <= ninputargs)

    flag = deblank(inputargs{narg});
    narg = narg + 1;
    % fprintf(1,'Argument: %s\n',flag);
    if(~isstr(flag))
      flag
      fprintf(1,'ERROR: non-string flag detected \n');
      s = []; return;
    end
    switch(flag)

      case '-h',
        if(argNcheck(flag,narg,ninputargs,1)) s = []; return; end
        s.hfig = inputargs{narg};
        narg = narg + 1;

      case '-vol',
        if(argNcheck(flag,narg,ninputargs,1)) s = []; return; end
        s.vol = inputargs{narg};
        if(length(size(s.vol)) < 3)
          fprintf('ERROR: vol has dimension less than 3\n');
          s = []; return;
        end
        s.redraw = 1;
        narg = narg + 1;

      case '-volid',
        if(argNcheck(flag,narg,ninputargs,1)) s = []; return; end
        s.volid = inputargs{narg};
        fprintf('INFO: loading %s\n',s.volid);
        s.vol = fmri_ldbvolume(s.volid);
        s.vol = permute(s.vol,[2 3 1 4]);
        s.redraw = 1;
        narg = narg + 1;

      case '-volres',
        if(argNcheck(flag,narg,ninputargs,1)) s = []; return; end
        s.volres = inputargs{narg};
        if(length(size(s.volres)) ~= 4)
          fprintf('ERROR: volres does not have 4 elements\n');
          s = []; return;
        end
        narg = narg + 1;

      case '-viewmosaic',
        s.mosaic = 1;

      case '-viewslice',
        s.mosaic = 0;

      case '-viewtoggle',
        s.mosaic = ~s.mosaic;

      case '-wbm',
        s.motion = 1;

      case '-kbd',
        s.curchar = get(s.hfig,'CurrentCharacter'); 
        fprintf('CurChar: %c\n',s.curchar);
        s = handlekbd(s,s.curchar);

      case '-wbd',
        s.mousedown = 1;
        xyz = get(gca,'CurrentPoint');
        c = round(xyz(1,1)) - 1;
        r = round(xyz(1,2)) - 1;
        newpos = s.curpos;
        newpos(1) = c;
        newpos(2) = r;
        if( inbounds(size(s.vol),newpos) )
          s.curpos = newpos;
	  v = s.vol(s.curpos(2)+1,s.curpos(1)+1,s.curpos(3)+1,s.curpos(4)+1);
          fprintf('%d %d %d %d  %g\n',s.curpos,v);
        end

      case '-setpos',
        if(argNcheck(flag,narg,ninputargs,1)) s = []; return; end
        newpos = inputargs{narg};
        narg = narg + 1;
        if( inbounds(size(s.vol),newpos) )
          if(s.curpos(3) ~= newpos(3) | s.curpos(4) ~= newpos(4)) 
            s.redraw = 1; end
          s.curpos = newpos;
	  v = s.vol(s.curpos(2)+1,s.curpos(1)+1,s.curpos(3)+1,s.curpos(4)+1);
          fprintf('%d %d %d %d  %g\n',s.curpos,v);
        end

      case '-wbu',
        s.mousedown = 0;

      otherwise
        fprintf(2,'ERROR: Flag %s unrecognized\n',flag);
        s = [];
        return;

    end % --- switch(flag) ----- %

  end % while(narg <= ninputargs)

return;
%--------------------------------------------------%

%--------------------------------------------------%
%% Check that there are at least N or more arguments %%
function r = argNcheck(flag,nflag,nmax,N)
  r = 0;
  if(nflag+N-1 > nmax) 
    fprintf(1,'ERROR: Flag %s needs %d argument(s)',flag,N);
    r = 1;
  end
return;

%--------------------------------------------------%
%% Check argument consistency, etc %%%
function s = check_params(s)

  if(isempty(s.hfig) & isempty(s.vol) & isempty(s.volid))
    fprintf('ERROR: must specify either vol or volid when init\n');
    s = []; return;
  end

  if(~isempty(s.vol) & ~isempty(s.volid))
    fprintf('ERROR: cannot specify both vol and volid \n');
    s = []; return;
  end

return;

%--------------------------------------------------%
%% Print Usage 
function print_usage
  fprintf(1,'USAGE: fast_volview\n');
  fprintf(1,' -h figurehandle \n');
  fprintf(1,' -vol vol4d (c,r,s,f)\n');
  fprintf(1,' -volid volid\n');
  fprintf(1,' -volres dc dr ds df \n');
  fprintf(1,' -title title \n');
  fprintf(1,' -setpos c r s f \n');
  fprintf(1,' -viewmosaic \n');
  fprintf(1,' -viewslice \n');
  fprintf(1,' -viewtoggle \n');

return
%--------------------------------------------------%

%--------------------------------------------------%
%% Print data structure
function s = isxavg_fe_print_struct(s,fid)
  if(nargin == 1) fid = 1; end

  ninvols = size(s.invols,1);
  fprintf(fid,'ninvols    %d\n',ninvols);
  for n = 1:ninvols
    if(~s.weighted)
      fprintf(fid,'  invol      %s\n',s.invols(n,:));
    else
      fprintf(fid,'  invol      %s  %g\n',s.invols(n,:),s.weights(n));
    end
  end
  fprintf(fid,'outvol      %s\n',s.avgvol);
  fprintf(fid,'pctsigch    %d\n',s.pctsigch);
  fprintf(fid,'firstslice  %d\n',s.firstslice);
  fprintf(fid,'nslices     %d\n',s.nslices);
  fprintf(fid,'logfile     %s\n',s.logfile);
  fprintf(fid,'synth       %s\n',s.synth);

return;
%--------------------------------------------------%

%--------------------------------------------------%
function s = initfigure
  s = main_struct;

  s.hfig   = figure;
  s.himage = imagesc(randn(64,64));
  s.haxis  = gca;
  set(s.himage,'EraseMode','none');
  set(gcf,'pointer','crosshair');
  set(gcf,'KeyPressFcn',          'fast_volview(''-kbd'');');
  set(gcf,'WindowButtonDownFcn',  'fast_volview(''-wbd'');');
  set(gcf,'WindowButtonUpFcn',    'fast_volview(''-wbu'');');
  %set(gcf,'WindowButtonMotionFcn','fast_volview(''-wbm'');');

  s.redraw = 1;
  set(gcf,'UserData',s);

return;

%--------------------------------------------------%
function r = inbounds(szvol,pos)
  r = 1;

  if(length(szvol) < 3) nslices = 1;
  else                  nslices = szvol(3);
  end

  if(length(szvol) < 4) nframes = 1;
  else                  nframes = szvol(4);
  end

  if(pos(1) < 0 | pos(1) >= szvol(2) | ...
     pos(2) < 0 | pos(2) >= szvol(1) | ...
     pos(3) < 0 | pos(3) >= nslices | ...
     pos(4) < 0 | pos(4) >= nframes) r = 0; end
return;

%--------------------------------------------------%
function s = handlekbd(s,c)

  switch(c)
    case {'+','='}
      s.curpos(3) = s.curpos(3) + 1;
      if(s.curpos(3) >= size(s.vol,3)) s.curpos(3) = 0; end
      s.redraw = 1;
    case {'-','_'}
      s.curpos(3) = s.curpos(3) - 1;
      if(s.curpos(3) < 0 ) s.curpos(3) = size(s.vol,3) - 1; end
      s.redraw = 1;
    case {'f'}
      s.curpos(4) = s.curpos(4) + 1;
      if(s.curpos(4) >= size(s.vol,4)) s.curpos(4) = 0; end
      s.redraw = 1;
    case {'b'}
      s.curpos(4) = s.curpos(4) - 1;
      if(s.curpos(4) < 0 ) s.curpos(4) = size(s.vol,4) - 1; end
      s.redraw = 1;
  end 

return;

