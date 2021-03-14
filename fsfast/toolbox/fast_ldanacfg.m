function cfg = fast_ldanacfg(fsfcfgfile)

% ldfsfcfg.m


%
% fast_ldanacfg.m
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

fsfcfgfile = 'fsf.cfg';
cfg = [];

fid = fopen(fsfcfgfile);
if(fid == -1) 
  fprintf('ERROR: could not open %s\n',fsfcfgfile);
  return;
end

% Scroll through blank/comment lines at beginning

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

cfg = fast_anacfgstruct;

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
   
   case 'usetpexclude'
     [cfg.usetpexclude n] = sscanf(line,'%*s %d',1);  
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
     [cfg.inorm n] = sscanf(line,'%*s %f',1);  
     if(n == 0) err = 1; end
    
   case 'runlistfile'
     [cfg.runlistfile n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
    
   case 'prewhiten'
     [cfg.prewhiten n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
    
   case 'evschrname'
     [cfg.evschrname n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
    
   case 'rfxextregrname'
     [rfxextregrname n] = sscanf(line,'%*s %s',1);  
     if(n ~= 0) 
       cfg.rfxextregrname = strvcat(cfg.rfxextregrname, rfxextregrname);
       % Read possible nextreg, not an error if not there %
       [nrfxextreg n] = sscanf(line,'%*s %*s %d',1);  
       nth = size(cfg.extregrname,1);
       if(n ~= 0) cfg.nrfxextreg(nth) = nrfxextreg; 
       else       cfg.nrfxextreg(nth) = -1;
       end
     else
       err = 1; 
     end
    
   case 'cond'
    ccfg = fast_readcondline(line);
    if(isempty(ccfg)) 
      err = 1; 
    else
      nthcond = length(cfg.cond) + 1;
      if(nthcond == 1) cfg.cond = ccfg;
      else             cfg.cond(nthcond) = ccfg;
      end
    end
    
   case 'maskstem'
     [cfg.maskstem n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
    
   case 'slicetiming'
     [cfg.slicetiming n] = sscanf(line,'%*s %s',1);  
     if(n == 0) err = 1; end
    
  end % switch
  
  if(err)
    fprintf('ERROR: tag %s, file %s\n',tag,cfgfile);
    fprintf('  LINE (%d):  %s\n',lineno,line);
    fclose(fid);
    cfg = [];
    return;
  end

end % while 

fclose(fid);

return;
%--------------------------------------------------------------%
