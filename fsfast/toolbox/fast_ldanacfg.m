function cfg = fast_ldanacfg(fsfcfgfile)

% ldfsfcfg.m

fsfcfgfile = 'fsf.cfg';
cfg = [];

fid = fopen(fsfcfgfile);
if(fid == -1) 
  fprintf('ERROR: could not open %s\n',fsfcfgfile);
  return;
end

% Scroll through blank/comment lines


% Read the first line %
line = fgets(fid);
if(line == -1) 
  fprintf('ERROR: %s is empty\n',fsfcfgfile);
  fclose(fid);
  return;
end
[tag n] = sscanf(line,'%s',1);  
if(n == 0) 
  fprintf('ERROR: %s is not formated properly\n',fsfcfgfile);
  fclose(fid);
  return;
end
if(~strcmpi(tag,'fsfastcfg'))
  fprintf('ERROR: %s is not formated properly\n',fsfcfgfile);
  fclose(fid);
  return;
end
[cfgversion n] = sscanf(line,'%*s %d',1);  
if(n == 0) 
  fprintf('ERROR: cannot find version number in %s \n',fsfcfgfile);
  fclose(fid);
  return;
end
if(cfgversion ~= 2) 
  fprintf('ERROR: version = %d, must be 2 (%s) \n',cfgversion,fsfcfgfile);
  fclose(fid);
  return;
end

cfg.version = 2;
cfg.erm = [];
cfg.extreg = '';

lineno = 0;
while(1)
  line = fgetl(fid);
  lineno = lineno + 1;
  if(isempty(line)) continue; end

  if(line == -1) break; end % End of file

  % Get the tag %
  [tag n] = sscanf(line,'%s',1);  

  if(n == 0) continue; end           % Empty line
  if(tag(1) == '#') continue; end    % Comment line

  %fprintf('tag = %s\n',tag);

  % Convert tag to lower case and switch %
  tag = lower(tag);
  err = 0;
  nexp = 1;
  switch(tag)
   
   case 'analysis'
     [cfg.analysis n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
   
   case 'polyfit'
     [cfg.polyfit n] = sscanf(line,'%*s %d',1);  
     if(n == 0) err = 1; end
    
   case 'nskip'
     [cfg.nskip n] = sscanf(line,'%*s %d',1);  
     if(n == 0) err = 1; end
   
   case 'useevschweight'
     [cfg.useevschweight n] = sscanf(line,'%*s %d',1);  
     if(n == 0) err = 1; end
   
   case 'tr'
     [cfg.TR n] = sscanf(line,'%*s %f',1);  
     if(n == 0) err = 1; end
    
   case 'fsd'
     [cfg.fsd n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
    
   case 'funcstem'
     [cfg.funcstem n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
    
   case 'inorm'
     [cfg.TR n] = sscanf(line,'%*s %f',1);  
     if(n == 0) err = 1; end
    
   case 'runlistfile'
     [cfg.runlistfile n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
    
   case 'prewhiten'
     [cfg.prewhiten n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
    
   case 'evschfname'
     [cfg.evschfname n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
    
   case 'extreg'
     [extreg n] = sscanf(line,'%*s %s',1);  
     if(n ~= 0) 
       cfg.extreg = strvcat(cfg.extreg,extreg);
       % Read possible nextreg, not an error if not there %
       [nextreg n] = sscanf(line,'%*s %*s %d',1);  
       nth = size(cfg.extreg,1);
       if(n ~= 0) cfg.nextreg(nth) = nextreg; 
       else       cfg.nextreg(nth) = -1;
       end
     else
       err = 1; 
     end
    
   case 'erm'
    ntherm = length(cfg.erm);
    [evid n]      = sscanf(line,'%*s %d',1);  
    if(n == 0) err = 1; end
    [label n]     = sscanf(line,'%*s %*d %s',1);  
    if(n == 0) err = 1; end
    [modelname n] = sscanf(line,'%*s %*d %*s %s',1);  
    if(n == 0) err = 1; end
    [bcw n]       = sscanf(line,'%*s %*d %*s %*s %f',1);  
    if(n == 0) err = 1; end
    [psdmin n]    = sscanf(line,'%*s %*d %*s %*s %*f %f',1);  
    if(n == 0) err = 1; end
    [dpsd n]      = sscanf(line,'%*s %*d %*s %*s %*f %*f %f',1);  
    if(n == 0) err = 1; end
    [psdmax n]    = sscanf(line,'%*s %*d %*s %*s %*f %*f %*f %f',1);  
    if(n == 0) err = 1; end
    params = [];
    fmt0 = '%*s %*d %*s %*s %*f %*f %*f %*f';
    while(1)
      fmt = sprintf('%s %%f',fmt0);
      [ptmp n] = sscanf(line,fmt,1);
      if(n==0) break; end
      params = [params ptmp];
      fmt0 = sprintf('%s %%*f',fmt0);
    end
    % Check nparams vs expected nparams
    if(~err) 
      cfg.erm(ntherm+1).evid = evid;
      cfg.erm(ntherm+1).label = label;
      cfg.erm(ntherm+1).model = modelname;
      cfg.erm(ntherm+1).bcw = bcw;
      cfg.erm(ntherm+1).psdmin = psdmin;
      cfg.erm(ntherm+1).dpsd = dpsd;
      cfg.erm(ntherm+1).psdmax = psdmax;
      cfg.erm(ntherm+1).params = params;
    end
    
   case 'mask'
     [cfg.mask n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
    
   case 'slicetiming'
     [cfg.slicetiming n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
    
   case 'periodic'
     [cfg.periodic.basis n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
     [cfg.periodic.freq n] = sscanf(line,'%*s %*s %g',1);  
     if(n == 0) err = 1; end
     [cfg.periodic.phase n] = sscanf(line,'%*s %*s %*s %g',1);  
     if(n == 0) err = 1; end
    
  end % switch
  
  if(err)
    tagerror(tag,fsfcfgfile,line,lineno,nexp,n);
    fclose(fid);
    cfg = [];
    return;
  end

end % while 

fclose(fid);

return;
%--------------------------------------------------------------%
function tagerror(tag,cfgfile,line,lineno,nexp,nact)
fprintf('ERROR: tag %s, file %s\n',tag,cfgfile);
fprintf('  LINE (%d):  %s\n',lineno,line);
fprintf('  Expected %d, found %d\n',nexp,nact);
return;
