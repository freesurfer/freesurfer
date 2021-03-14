function r = editmask(varargin)
% r = editmask(varargin)
% Edit a functional mask


%
% editmask.m
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

version = 'editmask.m @FS_VERSION@';
r = 1;

%% Print usage if there are no arguments %%
if(nargin == 0)
  print_usage;
  return;
end

hcurrentfig = get(0,'CurrentFigure');

%--------------- Initialize if needed -----------------%
init = 0;
flag = deblank(varargin{1});
if(~isempty(strmatch(flag,'-init')) | isempty(hcurrentfig))
  init = 1;
  s = parse_args(varargin);
  if(isempty(s)) return; end
  s = check_params(s);
  if(isempty(s)) return; end

  % Load the base %
  fprintf('Loading base %s\n',s.baseid);
  s.base = fmri_ldbvolume(s.baseid);
  s.base = s.base(:,:,:,s.baseframe+1);
  if(isempty(s.base))
    fprintf('ERROR: cannot load %s\n',s.baseid);
    return;
  end
  s.volsize = size(s.base);

  % Convert to mosaic %
  s.base = vol2mos(permute(s.base,[2 3 1]));
  s.mossize = size(s.base);

  % Load the input mask (if specificied) or segment the base %
  if(~isempty(s.inmaskid))
    s.mask = fmri_ldbvolume(s.inmaskid);
    if(isempty(s.mask))
      fprintf('ERROR: could not load input mask %s\n',s.inmaskid);
      return;
    end
    s.mask = s.mask(:,:,:,1);
    s.mask = vol2mos(permute(s.mask,[2 3 1]));
    if(s.invertinmask) s.mask = ~s.mask; end
  else
    s.mask = zeros(size(s.base));
    if(s.segbase)
      meanbase = mean(reshape1d(s.base));
      ind = find(s.base > s.segthresh*meanbase);
      s.mask(ind) = 1;
      fprintf('Segmenting base %g %d\n',meanbase,length(ind));
    end
    if(s.invertinmask) s.mask = ~s.mask; end
    s.saveneeded = 1;
  end

  if(s.heqbase)
     fprintf('Equalizing Base Intesity Values\n');
     [histbase valbase] = hist(reshape1d(s.base),100);
     pbase = histbase/sum(histbase);
     pdfbase = cumsum(pbase);
     ix = max(find(pdfbase < .98));
     x = valbase(ix);
     ind = find(s.base > x);
     s.base(ind) = x;
  end

  % Rescale Base to be between 1 and (s.ncmap-1) %
  minbase = min(reshape1d(s.base));
  maxbase = max(reshape1d(s.base));
  s.base = floor((s.ncmap-2)*(s.base-minbase)/(maxbase-minbase)) + 1;
  % Create the color map
  s.cmap = gray(s.ncmap); % Start with gray
  % Set the last entry to be color;
  s.cmap(s.ncmap,:) = [.75 .75 0]; % raises a figure window here (why?)

  s.displayimg = s.base .* (~s.mask) + s.ncmap*s.mask;

  % Image the display image %
  s.himage = image(s.displayimg);
  set(s.himage,'EraseMode','none');
  axis image;
  s.hfig = gcf;
  s.haxis = gca;
  colormap(s.cmap);
  set(gcf,'pointer','crosshair');

  % Set up the call-back functions %
  set(gcf,'KeyPressFcn',          'editmask(''kbd'');');
  set(gcf,'WindowButtonDownFcn',  'editmask(''wbd'');');
  set(gcf,'WindowButtonUpFcn',    'editmask(''wbu'');');
  set(gcf,'WindowButtonMotionFcn','editmask(''wbm'');');
  s.nmarked = length(find(s.mask==1));

  set(gcf,'UserData',s);

  fprintf('\n');
  fprintf(' ----------------------------------\n');
  fprintf(' For help press "h"\n');
  fprintf(' ------------ Help ----------------\n');
  printhelp;
  fprintf('\n');
  fprintf(' ------ Current State -------------\n');
  printstate(s);
  fprintf('\n');

  return;
end
%---------------------------------------------------------%
%---------------------------------------------------------%
%---------------------------------------------------------%

%----------- Parse the call-back function ----------%
s = get(gcf,'UserData');
redraw = 0;
switch(flag)

  case {'wbd','wbm'} % -------Window Button Down ------------ %
    if(strcmp(flag,'wbm') & ~s.mousedown) return; end
    if(s.editmode == 0) return; end

    xyz = get(gca,'CurrentPoint');
    c = round(xyz(1,1));
    r = round(xyz(1,2));
    if(r < 1 | r > size(s.base,1) | c < 1 | c > size(s.base,2))
      return;
    end
    s.cp = round([r c]);
    rlist = r-(s.brushsize-1):r+(s.brushsize-1);
    ind = find(rlist < s.mossize(1));
    rlist = rlist(ind);

    clist = c-(s.brushsize-1):c+(s.brushsize-1);
    ind = find(clist < s.mossize(2));
    clist = clist(ind);

    %fprintf('%d %d %g  %d\n',r,c,s.base(r,c),s.mask(r,c));

    for r = rlist,
      for c = clist,
        switch(s.editmode)
        case 1, s.mask(r,c) = 1;
        case 2, s.mask(r,c) = 0;
        case 3, s.mask(r,c) = ~s.mask(r,c);
        end
      end
    end

    s.displayimg = s.base .* (~s.mask);
    s.displayimg = s.base .* (~s.mask) + s.ncmap * s.mask;

    s.nmarked = length(find(s.mask==1));

    s.saveneeded = 1;
    s.mousedown = 1;
    redraw = 1; 

  case {'wbu'} % -------Window Button Up ------------ %
    %fprintf('Button Up\n');
    s.mousedown = 0;

  case {'kbd'} %----------- Keyboard -------------%
    c = get(s.hfig,'CurrentCharacter'); 
    switch(c)
      case {'b'},
        s.brushsize = s.brushsize + 1;
        if(s.brushsize > 8) s.brushsize = 1; end
        tmp = 2*(s.brushsize-1) + 1;
        fprintf('Brush Size: %d x %d\n',tmp,tmp);
        set(gcf,'UserData',s);
        return;
      case {'h'},
        fprintf('\n');
        printhelp;
        fprintf('\n');
        printstate(s);
        fprintf('\n');
      case {'m'},
        s.editmode = s.editmode + 1;
        if(s.editmode > 3) s.editmode = 0; end
        fprintf('Mode: %s\n',modename(s.editmode));
      case {'q'},
        if(s.saveneeded)
          resp = questdlg('Changes pending. Do you want to quit?',...
                          'Quit?','Quit','No','No');
          if(strcmp(resp,'No')) return; end
        end
        close(s.hfig);
        r = 0;
        return;
      case {'s'},
        if(isempty(s.outmaskid))
          curdir = pwd;
          cd(s.savedir);
          [fname pname] = uiputfile(s.outmaskid,'Save Mask');
          cd(curdir);
          if(fname ~= 0)
            if(isempty(pname)) ud.savedir  = '.';
            else               ud.savedir  = pname;
            end
            s.outmaskid = strcat(pname,fname);
          else return;  
          end
        end
        tmp = mos2vol(s.mask,[s.volsize(2) s.volsize(3) s.volsize(1)]);
        tmp = permute(tmp,[3 1 2]);
        fprintf('Saving to %s\n',s.outmaskid);
        fmri_svbvolume(tmp,s.outmaskid);
        fprintf('Done\n');
        clear tmp;
        s.saveneeded = 0;
        set(gcf,'UserData',s);
        return;
      case {'t'},
        s.showbaseonly = ~ s.showbaseonly;
        if(s.showbaseonly) 
          s.displayimg = s.base;
        else
          s.displayimg = s.base .* (~s.mask);
          s.displayimg = s.base .* (~s.mask) + s.ncmap*s.mask;
        end
        redraw = 1;
        figure(s.hfig);
        axes(s.haxis);
      case {'v'},
        fprintf('\n');
        fprintf('\n');
        printstate(s);
        fprintf('\n');
        fprintf('\n');
      case {'z'},
        zoom;
        s.zoomstate = ~ s.zoomstate;
        set(gcf,'UserData',s);
        if(s.zoomstate)
          fprintf('Zoom is on\n');
        else
          fprintf('Zoom is off\n');
        end
        return;
    end

end % --- switch(flag) ----- %

if(redraw)
  set(s.himage,'CData',s.displayimg);
end

set(gcf,'UserData',s);

return;
%--------------------------------------------------%
%--------------------------------------------------%
%--------------------------------------------------%

%--------------------------------------------------%
%% Print Usage 
function print_usage(dummy)

  fprintf(1,'USAGE:\n');
  fprintf(1,'  editmask\n');
  fprintf(1,'     -base     base volume  \n');
  fprintf(1,'     -heqbase  equalize base volume intensities \n');
  fprintf(1,'     -segbase  initial mask from segmented base volume  \n');
  fprintf(1,'     -segthresh use threshold to segment base volume  \n');
  fprintf(1,'     -inmask   load and edit this mask\n');
  fprintf(1,'     -outmask  set the output volume\n');
  fprintf(1,'     -invert (invert input mask or segbase)\n');
return
%--------------------------------------------------%

%--------------------------------------------------%
%% Default data structure
function s = editmask_struct
  s.volsize        = [];
  s.mossize        = [];
  s.baseid         = '';
  s.base           = [];
  s.heqbase            = 0;
  s.inmaskid       = '';
  s.invertinmask   = 0;
  s.outmaskid      = '';
  s.mask           = [];
  s.baseframe      = 0;
  s.hfig           = [];
  s.haxis          = [];
  s.himage         = [];
  s.ncmap          = 64;
  s.displayimg     = [];
  s.segbase        = 0;
  s.segthresh      = .75;
  s.cp             = [1 1]; % current point [r c]
  s.zoomstate      = 0;
  s.brushsize      = 2;
  s.mousedown      = 0;
  s.showbaseonly   = 0;
  s.savedir        = '.';
  s.editmode       = 0;
  s.nmarked        = 0;
  s.saveneeded     = 0;
return;

%--------------------------------------------------%
% ----------- Parse Input Arguments ---------------%
function s = parse_args(varargin)

  fprintf(1,'Parsing Arguments \n');
  s = editmask_struct;

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

      case {'-b','-base'},
        arg1check(flag,narg,ninputargs);
        s.baseid = inputargs{narg};
        narg = narg + 1;

      case {'-bf','-baseframe'}
        arg1check(flag,narg,ninputargs);
        s.baseframe = sscanf(inputargs{narg},'%d',1);
        narg = narg + 1;

      case {'-inmask'},
        arg1check(flag,narg,ninputargs);
        s.inmaskid = inputargs{narg};
        narg = narg + 1;

      case {'-invert'},
        s.invertinmask = 1;

      case {'-outmask'},
        arg1check(flag,narg,ninputargs);
        s.outmaskid = inputargs{narg};
        narg = narg + 1;

      case {'-segthresh'},
        arg1check(flag,narg,ninputargs);
        s.segthresh = sscanf(inputargs{narg},'%f');
        s.segbase = 1;
        narg = narg + 1;

      case '-segbase',
        s.segbase = 1;

      case '-heqbase',
        s.heqbase = 1;

      case '-verbose',
        s.verbose = 1;

      case {'-monly', 'umask'},% ignore
        arg1check(flag,narg,ninputargs);
        narg = narg + 1;

      case {'-debug','-echo','-init'}, % ignore

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
%% Check argument consistency, etc %%%
function s = check_params(s)
  %fprintf(1,'Checking Parameters\n');

  if( isempty(s.baseid) )
    fprintf(2,'ERROR: No base volumes specified\n');
    s=[]; return;
  end
return;

%--------------------------------------------------%
function name = modename(modeno)
  switch(modeno)
  case 0, name = 'No Edit'; 
  case 1, name = 'Set'; 
  case 2, name = 'Unset'; 
  case 3, name = 'Toggle'; 
  end
return

%--------------------------------------------------%
function printstate(s)
  s.nmarked = length(find(s.mask==1));
  tmp = 2*(s.brushsize-1) + 1;
  fprintf('Base Volume: %s\n',s.baseid);
  fprintf('Input  Mask Volume: %s\n',s.inmaskid);
  fprintf('Output Mask Volume: %s\n',s.outmaskid);
  fprintf('Total number of voxels:  %d\n',prod(s.volsize));
  fprintf('Number of voxels marked: %d\n',s.nmarked);
  fprintf('Brush Size: %d x %d\n',tmp,tmp);
  fprintf('Edit Mode: %s\n',modename(s.editmode));
  if(s.saveneeded)
    fprintf('Changes Pending\n');
  else
    fprintf('No Changes Pending\n');
  end
return;

%--------------------------------------------------%
function printhelp
  fprintf('b - change brush size (1x1,3x3,5x5,7x7,9x9,11x11,13x13,...)\n');
  fprintf('h - print help (this message)\n');
  fprintf('m - change edit mode\n');
  fprintf('    No Edit - button down does nothing\n');
  fprintf('    Set     - button down turns on mask at voxel\n');
  fprintf('    Unset   - button down turns off mask at voxel\n');
  fprintf('    Toggle  - button down toggels mask at voxel\n');
  fprintf('q - quit/exit\n');
  fprintf('s - save mask\n');
  fprintf('t - toggle mask on and off\n');
  fprintf('v - print the current state\n');
  fprintf('z - toggle zoom state. When zoom is on, use left\n');
  fprintf('    button to zoom in and right to zoom out\n');
return;
