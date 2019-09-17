pad = 1 ;
[ivol,M,mr] = load_mgh('EC4_40x_25um_crop.Deformed.mgz') ;
ivol = ivol - mean(ivol(:)) ;

disp(sprintf('computing neural convolutional layer')) ;
[neuron_vol,M2,mr2] = load_mgh('EC4.neuron_label2.mgz') ;
ind = find(neuron_vol > 0) ;
[subx, suby] = ind2sub(size(neuron_vol), ind);
xmin = min(subx) ; xmax = max(subx) ;
ymin = min(suby) ; ymax = max(suby) ;
template_neuron = ivol(xmin-pad:xmax+pad,ymin-pad:ymax+pad) ;
template_neuron = template_neuron ./ norm(template_neuron) ;
xc_neuron = conv2(ivol, template_neuron,'same');

disp(sprintf('computing fiber convolutional layer')) ;
[fiber_vol,M2,mr2] = load_mgh('EC4.fiber_label.mgz') ;
ind = find(fiber_vol > 0) ;
[subx, suby] = ind2sub(size(fiber_vol), ind);
xmin = min(subx) ; xmax = max(subx) ;
ymin = min(suby) ; ymax = max(suby) ;
template_fiber = ivol(xmin-pad:xmax+pad,ymin-pad:ymax+pad) ;
template_fiber = template_fiber ./ norm(template_fiber) ;
xc_fiber = conv2(ivol, template_fiber,'same')/sum(template_fiber(:));


disp(sprintf('computing glial convolutional layer')) ;
[glia_vol,M2,mr2] = load_mgh('EC4.glia_label.small.mgz') ;
ind = find(glia_vol > 0) ;
[subx, suby] = ind2sub(size(glia_vol), ind);
xmin = min(subx) ; xmax = max(subx) ;
ymin = min(suby) ; ymax = max(suby) ;
template_glia = ivol(xmin-pad:xmax+pad,ymin-pad:ymax+pad) ;
template_glia = template_glia / norm(template_glia);
xc_glia = conv2(ivol, template_glia,'same');

disp(sprintf('computing background convolutional layer')) ;
[bg_vol,M2,mr2] = load_mgh('EC4.bg_label.mgz') ;
ind = find(bg_vol > 0) ;
[subx, suby] = ind2sub(size(bg_vol), ind);
xmin = min(subx) ; xmax = max(subx) ;
ymin = min(suby) ; ymax = max(suby) ;
template_bg = ivol(xmin-pad:xmax+pad,ymin-pad:ymax+pad) ;
template_bg = template_bg ./ norm(template_bg) ;
xc_bg = conv2(ivol, template_bg,'same');

%xc_neuron = xc_neuron - (xc_fiber + xc_glia) ;
%xc_fiber = xc_fiber - (xc_neuron + xc_glia) ;
%xc_glia = xc_glia - (xc_fiber + xc_neuron) ;

disp(sprintf('saving neurons')) ;
save_mgh(xc_neuron, 'xcn.mgz', M,mr) ;
disp(sprintf('saving fibers')) ;
save_mgh(xc_fiber, 'xcf.mgz', M,mr) ;
disp(sprintf('saving glia')) ;
save_mgh(xc_glia, 'xcg.small.mgz', M,mr) ;
disp(sprintf('saving background')) ;
save_mgh(xc_bg, 'xcb.mgz', M,mr) ;




if 0
%xc = conv2(ivol, rot90(conj(template),2),'same');
xc = filter2(template, ivol) ;

xcfull = xcorr2(template, ivol) ;
[trows,tcols] = size(template) ;
[irows,icols] = size(neuron_vol) ;
first_row = floor(trows/2) ;
last_row = first_row+irows-1 ;
first_col = floor(tcols/2) ;
last_col = first_col+icols-1 ;
xc = xcfull(first_row:last_row,first_col:last_col) ;
end
