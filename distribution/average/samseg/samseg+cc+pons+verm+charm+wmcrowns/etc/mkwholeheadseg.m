% 

% GTM Seg Pons belly, vermis

% Pituitary
% Sclimbic: Fornix AC

% erase 258 (other tissue)
% erase 165 (skull)

% Charm
% 501       Air-Internal  -->     262 Sinus (merge 0)
% 502       Artery  --> 502 (merge 0)  
% 506       EyeBalls  259 (265,266)
% 507       Other-Tissues   - keep (merge with CSF?)
% 508       Rectus-Muscles         
% 509       Mucosa (merge 0)
% 511       Skin   (merge 0?)
% 512       Spinal-Cord    --> 126
% 514       Vein     merge 24=CSF 8,47=clbumctx      
% 515       Bone-Cortical   keep zero:165 merge 24=CSF
% 516       Bone-Cancellous keep zero:165 merge 24=CSF
%517       Background          
% 520       Cortical-CSF       -->24
% 530       Optic-Nerve        85=chiasm?

cd /autofs/cluster/fssubjects/atlases/fsv6.aseg_atlas
slist = char(textread('fsv6.aseg.train.slist.txt','%s'));
%slist = 'v6.vc763';
nslist = size(slist,1);

% Center of the cblum (for vermis)
mni305cblum_crs = [127.5, 158, 65];

tic;
for n = 1:nslist
  subject = deblank(slist(n,:));
  fprintf('%d %s %6.1f ===================================\n',n,subject,toc);

  fname = sprintf('%s/mri/seg_edited10.mgz',subject);
  manseg = MRIread(fname);
  fname = sprintf('%s/mri/aseg.mgz',subject);
  aseg = MRIread(fname);
  fname = sprintf('%s/mri/aparc.a2009s+aseg.mgz',subject);
  a2009s = MRIread(fname);
  fname = sprintf('%s/mri/apas+head.mgz',subject);
  ponsverm = MRIread(fname);
  fname = sprintf('%s/mri/charm.bcd/seg.mgz',subject);
  charm = MRIread(fname);
  fname = sprintf('%s/mri/samseg.bcd/seg.mgz',subject);
  samseg = MRIread(fname);
  fname = sprintf('%s/mri/seg_edited10+wmcrowns.mgz',subject);
  wmcrowns = MRIread(fname); % lh=34, rh=66

  
  fname = sprintf('%s/mri/transforms/talairach.xfm.lta',subject);
  [subj2mni305 lta] = lta_read(fname); % vox2vox
  center_crs0 = inv(subj2mni305)*[127.5 127.5 127.5 1]';% mnicenter
  center_col = center_crs0(1) + 1;
  cblum_crs0 = inv(subj2mni305)*[mni305cblum_crs 1]'; % cblum center
  cblum_center_col = cblum_crs0(1) + 1;

  % Not all subjects have optic chiasm or left choroid. In these
  % cases, get the segs from the aseg, just so we can have
  % data set where all subjects have all ROIs
  ind = find(manseg.vol == 85);  % Chiasm
  if(length(ind)==0)
    ind = find(aseg.vol == 85);  
    manseg.vol(ind) = 85;
  end
  ind = find(manseg.vol == 31);  
  if(length(ind)==0)
    ind = find(aseg.vol == 31);  
    manseg.vol(ind) = 31;
  end
  
  % handle "undetermined". I can only find a left undet in one
  % subject v6.991102_vc1401.  It appears to be some kind of
  % stroke, so just say it it WMSA
  ind = find(manseg.vol == 29);
  if(length(ind)>0) manseg.vol(ind)=77; end
  ind = find(manseg.vol == 61);
  if(length(ind)>0) manseg.vol(ind)=77; end
  
  % handle "5th vent". This is only in 11 of the aseg subjects and
  % never exceeds 30 voxels. It is bet the LH and RH ventricles at
  % the most anterior end. It does not appear to be anything in the
  % subjects I looked at.
  ind = find(manseg.vol == 72);
  if(length(ind)>0) 
    [rr cc ss] = ind2sub(manseg.volsize,ind);
    ind2 = find(cc<=center_col);
    manseg.vol(ind(ind2)) =  4; % left lat vent
    ind2 = find(cc>center_col);
    manseg.vol(ind(ind2)) = 43; % right lat vent
  end
  
  % Vessel 30 62. All cases have, but some are very small. These
  % are likely VR spaces. Since all cases have both left and right, 
  % keep these in the atlas.

  % non-WM-hypointensities 80 (81,82); 12 cases have not systematic,
  % some in putamen, some in caud, some in cblum, distributed even in
  % one subject. To fix will need to create clusters, then look around
  % each cluster and fill in. What a pain.
  
  newseg = manseg;
  
  % CC - set to 192 to reduce the number of segments to speed up.
  % In FS, the individual segs are just generated using a PCA.
  for id = 251:255
    ind = find(aseg.vol==id);
    newseg.vol(ind) = 192; % id;
  end
  
  % Pons 
  id = 174;
  ind = find(ponsverm.vol==id);
  newseg.vol(ind) = id;
  
  % Vermis - separate into left and right
  % 6/18/2023 - unfortunately, just noticed that I got left and
  % right backwards. Will reverse in the LUT so I don't have
  % to regenerate volumes and atlas
  id = 172
  ind = find(ponsverm.vol==id);
  [rr cc ss] = ind2sub(ponsverm.volsize,ind);
  ind2 = find(cc<=cblum_center_col);
  newseg.vol(ind(ind2)) = 183; % left verm
  ind2 = find(cc>cblum_center_col);
  newseg.vol(ind(ind2)) = 184; % right verm

  % High Myelin (L=11k, R=12k) - these match the apas in v6
  % 129  ctx-lh-G_precentral
  % 146  ctx-lh-S_central
  % 128  ctx-lh-G_postcentral
  % 133  ctx-lh-G_temp_sup-G_temp_transv_and_interm_S
  % 143  ctx-lh-Pole_occipital
  % 111  ctx-lh-G_cuneus
  % 145  ctx-lh-S_calcarine
  % 122  ctx_lh_G_oc-temp_med-Lingual
  for hemi = 1:2
    ctxhimyid = 11300 + 1000*(hemi-1);
    ntot = 0;
    for roi = [129 146 128 133 143 111 145 122]
      id = 11000 + roi + 1000*(hemi-1);
      ind = find(a2009s.vol == id);
      newseg.vol(ind) = ctxhimyid;
      ntot = ntot + length(ind);
    end
    fprintf('id %d n=%d\n',ctxhimyid,ntot);
  end
  
  % WM Crowns
  for id = [34 66]
    ind = find(wmcrowns.vol==id);
    newseg.vol(ind) = id;
  end
  
  % Charm
  idlist = [501 502 506 507 508 509 511 512 514 515 516 520 530];
  for id = idlist
    ind = find(charm.vol==id & ~manseg.vol);
    if(id==501)     newseg.vol(ind) = 262; % sinus
    elseif(id==506) newseg.vol(ind) = 259; % eye
    %elseif(id==507) newseg.vol(ind) = 258; % head-extracerebral tricky
    %elseif(id==515) newseg.vol(ind) = 165; % skull
    elseif(id==520) newseg.vol(ind) = 24;  % csf
    elseif(id==512) newseg.vol(ind) = 126;  % spinal cord
    else newseg.vol(ind) = id+400;
    end
    ulist = unique(samseg.vol(ind));
    fprintf('id %4d ',id)
    fprintf('%4d ',ulist)
    fprintf('\n');
  end
  
  % Replace voxels that are cerebral or cerebellar cortex in charm but
  % 0 in manseg with CSF
  ind=find((charm.vol==3 | charm.vol==42 | ...
	    charm.vol==8 | charm.vol==47) & ~manseg.vol);
  newseg.vol(ind) = 24;
  
  % Extend brainstem if needed
  ind=find(charm.vol==16 & ~manseg.vol);
  newseg.vol(ind) = 16;
  
  fname = sprintf('%s/mri/seg_edited10+cc+pons+verm+charm+wmcrowns.22.12.31.mgz',subject);
  MRIwrite(newseg,fname);

end

