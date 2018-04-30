use_fisher = 1 ;
use_air = 0 ;

str = getenv('fluidthresh');
if (length(str) > 0)
  fluidthresh = sscanf(str, '%f')
else
  fluidthresh = 30
end

base = getenv('base');

recon=sprintf('%s',base);
mdir=sprintf('%s/mri',recon) ;
sdir=sprintf('%s/scripts',recon)' ;
ldir=sprintf('%s/label',recon) ;
odir=sprintf('%s/mri/orig',base) ;
datdir=sprintf('%s/mri/orig/opt_files',base) ;


wm = 1 ;
gm = 2 ;
fluid = 3 ;
air = 4 ;
labels = str2mat('wm', 'gm', 'fluid', 'air') ;
labels = str2mat('wm', 'gm', 'fluid') ;
nlabels = size(labels,1) ;
T1 = 1 ; PD =2 ; T2star = 3 ;
vol_names = str2mat('T1', 'PD', 'T2star') ;
nvols = size(vol_names,1) ;


disp(sprintf('loading labeled data...')) ;

flist = sprintf('%s/flist.dat', datdir) ;;
fid = fopen(flist, 'r') ;
if (fid < 0)
  error(sprintf('could not open file list %s', flist));
end


n = 1 ;
line = fgetl(fid) ;
while (line >= 0)
  slash_ind = strfind(line, '/');
  if ( length(slash_ind > 0))
    fname = line(slash_ind(end)+1:end);
  else
    fname = line;
  end
  disp(sprintf('%d: loading label vals for %s', n,line));
  for lno=1:nlabels
      vals = load(sprintf('%s/%s.%s.dat', datdir, fname, deblank(labels(lno,:))));
      switch (lno)
           case wm
            wm_vals(n,:) = vals ;
           case gm
            gm_vals(n,:) = vals ;
           case air
            air_vals(n,:) = vals ;
           case fluid
            fluid_vals(n,:) = vals ;
          end
  end
  n = n + 1 ;
  line = fgetl(fid) ;
end
fclose(fid) ;

%normalize the observations
nwm_vals = wm_vals ;
ngm_vals = gm_vals ;
if (use_air)
  nair_vals = air_vals ;
end
nfluid_vals = fluid_vals ;
for r=1:size(wm_vals,1)
    gm_norms(r) = norm(gm_vals(r,:)) ;
    wm_norms(r) = norm(wm_vals(r,:)) ;
    fluid_norms(r) = norm(fluid_vals(r,:)) ;
    nwm_vals(r,:) = wm_vals(r,:)./wm_norms(r) ;
    ngm_vals(r,:) = gm_vals(r,:)./gm_norms(r) ;
    nfluid_vals(r,:) = fluid_vals(r,:)./fluid_norms(r) ;
    if  (use_air)
      air_norms(r) = norm(air_vals(r,:)) ;
      nair_vals(r,:) = air_vals(r,:)./ air_norms(r) ;
    end
end

% normalize observations to prevent ill-conditioning do to differential scaling
% normalizing does help the condition number, but results in bright
% stuff just outside brain, so it is off for now
normalize = 0 ;
if (normalize)   
  mg = mean(ngm_vals') ;
  cg = cov(ngm_vals') ;
  mw = mean(nwm_vals') ;
  cw = cov(nwm_vals') ;
else
  wm_norms = ones(size(wm_norms)) ;
  gm_norms = ones(size(gm_norms)) ;
  mg = mean(gm_vals') ;
  cg = cov(gm_vals') ;
  mw = mean(wm_vals') ;
  cw = cov(wm_vals') ;
end

cinv = cg+cw;

%keyboard
while (cond(cinv) > 1e10)
  cinv = cinv + eye(size(cinv))*.01*mean(diag(cinv));
end

c = inv(cinv);

wgw = c * ((mw./wm_norms)' - (mg./gm_norms)') ;
wgw = c * ((mw)' - (mg)') ;
%wgw =  (mw' - mg') ;
%wgw = wgw ./ norm(wgw) ;

if (mean(wgw'*(mw'-mg')) < 0)
   wgw = wgw * -1 ;
end
offset =  mw*c*mw'-mg*c*mg';

% use nonlinear optimization to find weighting instead to reduce
% impact of fluid (and maybe air)
  if ((use_fisher == 0) & use_air)
  wgw = compute_optimal_weighting(wm_vals, gm_vals, fluid_vals, air_vals) ;
end


dark_vals = [fluid_vals] ;
bright_vals = [gm_vals, wm_vals] ;
md = mean(dark_vals') ;
mb = mean(bright_vals') ;
cd = cov(dark_vals') ;
cb = cov(bright_vals') ;

cinv = cb+cd ;
while (cond(cinv) > 1e10)
  cinv = cinv + eye(size(cinv))*.01*mean(diag(cinv));
end
c = inv(cinv) ;


wfmask = c * (mb' - md') ;
wfmask = wfmask ./ norm(wfmask) ;
if (mean(wfmask'*(mb'-md')) < 0)
   wfmask = wfmask * -1 ;
end

if  (use_air)
  dark_vals = [air_vals] ;
  bright_vals = [gm_vals, wm_vals] ;
  md = mean(dark_vals') ;
  mb = mean(bright_vals') ;
  cd = cov(dark_vals') ;
  cb = cov(bright_vals') ;

  cinv = cb+cd ;
  while (cond(cinv) > 1e10)
    cinv = cinv + eye(size(cinv))*.01*mean(diag(cinv));
  end
  c = inv(cinv) ;
  wamask = c * (mb' - md') ;
  wamask = wamask ./ norm(wamask) ;
  if (mean(wamask'*(mb'-md')) < 0)
     wamask = wamask * -1 ;
  end
end


clear('opt_vol') ;
n = 1 ;
fid = fopen(flist, 'r') ;
line = fgetl(fid) ;
while (line >= 0)
    if (exist('mri') ~= 1)   % read in vol with MRIread so we have the header
       mri = MRIread(sprintf('%s.mgz', line));
    end
    [v,M,mr]  = load_mgh(sprintf('%s.mgz', line)) ;
    if (exist('opt_vol') == 0)
       opt_vol = zeros(size(v)) ;
       fluid_mask_vol = zeros(size(v)) ;
        if (use_air)
         air_mask_vol = zeros(size(v)) ;
        end
    end
    opt_vol = opt_vol + wgw(n)*v ;
    fluid_mask_vol = fluid_mask_vol + wfmask(n)*v ;
    if (use_air)
      air_mask_vol = air_mask_vol + wamask(n)*v ;
    end
    n = n+1 ;
    line = fgetl(fid) ;
end
fclose(fid) ;


% using the Fisher discriminant will just optimize gray/white
% contrast, so  fluid and air may wind up bright and we need to mask
% them out
if (use_fisher)
%  mind = find(fluid_mask_vol <= 40) ;
%  opt_vol(mind) = zeros(size(mind));
  if (use_air)
    ind = find(air_mask_vol <= -50) ;
    opt_vol(ind) = zeros(size(ind));
    ind = find(air_mask_vol <= 100 & fluid_mask_vol < 0) ;
    opt_vol(ind) = zeros(size(ind));
  end
end

% keep track of zero locations and reset them to 0 later after we scale
%zind = find(opt_vol == 0);  

[gm_label,coords] = read_label('', sprintf('%s/gm.label',ldir));
ras = gm_label(:,2:5)';
ras(4,:) = ones(size(ras(4,:)));
if (strcmp(coords, 'scanner'))
  Mr2v = inv(mri.vox2ras);   % for scanner ras labels (new as of 11/1/2017)
else
  Mr2v = inv(mri.tkrvox2ras);   % for tkreg ras labels (old)
end

vox = Mr2v*ras;
gsub =  [vox(2,:)+1; vox(1,:)+1; vox(3,:)+1];
gsub =  round([vox(1,:)+1; vox(2,:)+1; vox(3,:)+1]);
gind = sub2ind(size(v), gsub(1,:), gsub(2,:),gsub(3,:)) ;

[wm_label,coords] = read_label('', sprintf('%s/wm.label', ldir));
ras = wm_label(:,2:5)';
ras(4,:) = ones(size(ras(4,:)));
if (strcmp(coords, 'scanner'))
  Mr2v = inv(mri.vox2ras);
else
  Mr2v = inv(mri.tkrvox2ras);
end

vox = Mr2v*ras;
wsub =  [vox(2,:)+1; vox(1,:)+1; vox(3,:)+1];
wsub =  round([vox(1,:)+1; vox(2,:)+1; vox(3,:)+1]);
wind = sub2ind(size(v), wsub(1,:), wsub(2,:),wsub(3,:)) ;

[fluid_label,coords] = read_label('', sprintf('%s/fluid.label',ldir));
ras = fluid_label(:,2:5)';
ras(4,:) = ones(size(ras(4,:)));
if (strcmp(coords, 'scanner'))
  Mr2v = inv(mri.vox2ras);
else
  Mr2v = inv(mri.tkrvox2ras);
end
vox = Mr2v*ras;
flsub =  [vox(2,:)+1; vox(1,:)+1; vox(3,:)+1];
flsub =  round([vox(1,:)+1; vox(2,:)+1; vox(3,:)+1]);
flind = sub2ind(size(v), flsub(1,:), flsub(2,:),flsub(3,:)) ;

mean_wm = mean(opt_vol(wind));
mean_gm = mean(opt_vol(gind));


% compute linear scaling that maps wm to 110 and gm to 50
%keyboard ;
m = 60 / (mean_wm - mean_gm) ;
b = 50- mean_gm * (60 / (mean_wm-mean_gm));
opt_vol = opt_vol * m + b ;
%opt_vol(zind) = zeros(size(zind)) ;
nind = find(opt_vol < 0 | opt_vol > 255) ;
nind = find(opt_vol < 0) ;
opt_vol(nind) = zeros(size(nind)) ;
opt_vol_unmasked = opt_vol ;

gm_mean = mean(fluid_mask_vol(gind)) ;
fl_mean = mean(fluid_mask_vol(flind)) ;
m = 100 / (gm_mean - fl_mean) ;
b = 100 - gm_mean * (100/(gm_mean-fl_mean));
fluid_mask_vol = fluid_mask_vol * m + b ;
nind = find(fluid_mask_vol < fluidthresh) ;
opt_vol(nind) = zeros(size(nind)) ;

Nopen = 1;
mind = find(fluid_mask_vol > fluidthresh);
nind = find(fluid_mask_vol <= fluidthresh);
mask = fluid_mask_vol;
mask(mind) = ones(size(mind));
mask(nind) = zeros(size(nind));
mask_opened = dilate3d(erode3d(mask,Nopen),Nopen);
fmask = fluid_mask_vol;                            
ind = find(mask_opened ==0);
fmask(ind) = zeros(size(ind));
ind = find(fmask < fluidthresh) ;
fmask(ind) = zeros(size(ind));

% don't let things get too bright and compress the rest of the values
hi_thresh = 300 ;
ind = find(fmask > hi_thresh) ;
fmask(ind) = hi_thresh*ones(size(ind));
ind = find(fmask == 0) ;
opt_vol(ind) = zeros(size(ind)) ;


disp(sprintf('saving  output volumes\n'));
save_mgh(fmask, sprintf('%s/brainmask.mgz', odir), M,mr);
save_mgh(opt_vol, sprintf('%s/opt.mgz', odir), M,mr);
save_mgh(opt_vol_unmasked, sprintf('%s/opt.unmasked.mgz', odir), M,mr);
save_mgh(fluid_mask_vol, sprintf('%s/fluid_mask.mgz', odir), M,mr);
if (use_air)
  save_mgh(air_mask_vol, sprintf('%s/air_mask.mgz', odir), M,mr);
end

if 0
PD_thresh = 1500 ;
T1_thresh = 750;


    wgw = compute_optimal_weighting(wm_vals, gm_vals, fluid_vals, air_vals) ;
    wgw = compute_optimal_pairwise_weighting(wm_vals, gm_vals) ;
    
    dark_vals = [gm_vals,fluid_vals,air_vals] ;
    out_vals = [fluid_vals,air_vals] ;
    wout = compute_optimal_pairwise_weighting(gm_vals, air_vals) ;
    wfluid = compute_optimal_pairwise_weighting(gm_vals, fluid_vals) ;
    disp(sprintf('loading volumes')) ;
    n = 1 ;
    w = wgw ./ norm(wgw) ;
    wo = wout ./ norm(wout) ;
    wf = wfluid ./ norm(wfluid) ;
    vol = zeros(size(v)) ;
    mask_vol = zeros(size(v)) ;
    fluid_mask = zeros(size(v)) ;
    if use_mef
        for f=flips
            disp(sprintf('loading echoes for flip %d', f));
            for e=echos
                [v,M,mr] = load_mgh(sprintf('mef%d_echo%d_avg.mgz', f, e)) ;
                %if you don't have multiple runs for each flip, and therefore don't use
                %the avg.csh script, run this line instead:
                %[v,M,mr] = load_mgh(sprintf('mef%d_echo%d.mgz', f, e)) ;
                vol = vol + w(n) * v ;
                mask_vol = mask_vol + wo(n) * v ;
                fluid_mask = fluid_mask + wf(n) * v ;
                n = n + 1 ;
            end
        end
    end
    
    
    for vn=1:nvols
        vol_name = remove_spaces(vol_names(vn,:)) ;
        [v,M,mr] = load_mgh(sprintf('../parameter_maps/%s.mgz', vol_name)) ;
        vol = vol + w(n) * v ;
        mask_vol = mask_vol + wo(n) * v ;
        fluid_mask = fluid_mask + wf(n) * v ;
        n = n+1 ;
        if (vn == T1)
            T1_vol = v ;
        elseif (vn == PD)
            PD_vol = v ;
        elseif (vn == T2star)
            T2star_vol = v ;
        end
    end
    vol = vol - min(vol(:)) ;
    vol_unmasked = vol ;
    mask_vol = mask_vol - min(mask_vol(:)) ;
    fluid_mask = fluid_mask - min(fluid_mask(:)) ;
    save_mgh(mask_vol, 'mask.mgz', M, mr);
    save_mgh(fluid_mask, 'fluid_mask.mgz', M, mr);
    
    % mask == 0 means that the voxel should be removed
    mask = ones(size(vol)) ;
    ind = find(T1_vol >T1_thresh) ;     % voxels with a T1 too long to be brain
    mask(ind) = zeros(size(mask(ind))) ;
    
    ind = find(PD_vol <PD_thresh) ;    % voxels with a PD too low to be brain
    mask(ind) = zeros(size(mask(ind))) ;
    
    mask_dilated = dilate3d(mask) ;  % restore border voxels
    mask_eroded = erode3d(mask) ;    % for finding border voxels in original mask
    
    % find voxels that are in the edge of the dilated mask but with really low PD and remove them
    ind_to_remove = find(mask_dilated > 0 & mask == 0 & PD_vol < (.75*PD_thresh));
    mask_dilated(ind_to_remove) = zeros(size(mask(ind_to_remove)));
    ind_to_remove = find(mask_eroded == 0 & mask ~= 0 & PD_vol <1.1*PD_thresh & T1_vol>350 & T2star_vol>50);
    mask_dilated(ind_to_remove) = zeros(size(mask(ind_to_remove)));
    ind_to_remove = find(mask_dilated == 1 & mask == 0 & PD_vol <1.1*PD_thresh & T1_vol>350 & T2star_vol>50);
    mask_dilated(ind_to_remove) = zeros(size(mask(ind_to_remove)));
    
    
    vol = vol_unmasked ;
    ind = find(mask_dilated == 0) ;
    vol(ind) = zeros(size(vol(ind))) ;
    save_mgh(mask, 'mask.mgz', M,mr) ;
    ind = find(vol > 0) ;
    mn = mean(vol(ind)) ;
    %vol = vol * 110/mn;
    save_mgh(vol, 'opt_gw.mgz', M, mr);
    
    unix(sprintf('mri_label_vals opt_gw.mgz %s/gm.label > opt_files/opt_gm.dat',ldir)) ;
    unix(sprintf('mri_label_vals opt_gw.mgz %s/wm.label > opt_files/opt_wm.dat',ldir)) ;
    
    gv = load('opt_files/opt_gm.dat') ;
    wv = load('opt_files/opt_wm.dat') ;
    x1 = mean(gv) ;
    x2 = mean(wv) ;
    y1 = 60 ;   % gm target
    y2 = 110 ;  % wm target
    m = (y2 - y1) / (x2 - x1) ;
    b = y1 - m * x1 ;
    svol = vol .* m + b ;
    ind = find(svol < 0) ;
    svol(ind) = zeros(size(ind)) ;
    ind = find(svol > 255) ;
    svol(ind) = 255*ones(size(ind)) ;
    
    mask = erode3d(mask_dilated) ;
    %save_mgh(svol, 'opt_gw.mgz', M, mr) ;
    ind_to_remove = find(mask_dilated > 0 & mask == 0 & T2star_vol>55 &svol >60 & PD_vol < 1.1*PD_thresh & T1_vol > 425);
    svol(ind_to_remove) = zeros(size(ind_to_remove)) ;
    ind_to_remove = find(mask_dilated > 0 & mask == 0 & T2star_vol>100 &svol >60  & T1_vol > 425);
    svol(ind_to_remove) = zeros(size(ind_to_remove)) ;
    ind_to_remove = find(mask_dilated > 0 & mask == 0 & PD_vol < PD_thresh &svol >60  & T1_vol > 425);
    svol(ind_to_remove) = zeros(size(ind_to_remove)) ;
    
    ind_to_remove = find(mask_dilated > 0 & mask == 0 & PD_vol <1.2*PD_thresh &svol >60  & T1_vol > 425 & T2star_vol <20);
    svol(ind_to_remove) = zeros(size(ind_to_remove)) ;
    
    save_mgh(svol, 'opt_gw.mgz', M, mr) ;
    
end
