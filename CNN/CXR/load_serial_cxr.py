import glob
import imageio
import pydicom
import gdcm
import os
import numpy as np
import surfa as sf
import pdb as gdb
import neuron as ne
from scipy import ndimage
from scipy.ndimage.interpolation import zoom
import scipy.ndimage.morphology as morph


def load_serial_cxr(base_path, subject_names='CVM*', study_names='CVA*'):
    subject_name_list = glob.glob(os.path.join(base_path, subject_names))
    nsubjects = len(subject_name_list)
    subject_list = []
    image_list = []
    for subject in subject_name_list:
        studies = glob.glob(os.path.join(subject, study_names))
        study_list = []
        imlist = []
        for study in studies:
            scans = glob.glob(os.path.join(study, '*'))

            # should only be 1 scan/study. Find the right one
            scan_list = []
            for scan in scans:
                dnames = glob.glob(os.path.join(scan,'*.dcm'))
                if len(dnames) < 1:
                    continue
                slist = []
                # more than 1 dicom in a dir means it was edge-enhanced 
                # or something went wrong - use the last one
                for dname in dnames:
                    im = pydicom.dcmread(dname)
                    if 'SECONDARY' in im.ImageType:
                        im = None
                        continue
                    
                    if im.SeriesDescription.find('AP')<0:
                        im = None
                        continue
                    if hasattr(im, 'DerivationDescription'):
                        if im.DerivationDescription.find('CATH') >= 0:
                            im = None
                            continue
                    slist.append(im)

                    im = slist[-1]
                if im is not None:
                    scan_list.append(im)
            if len(scan_list) > 0:
                im = scan_list[-1]

            # could check im.SeriesTime to pick last one
            if im is not None:
                imlist.append(im)
                study_list.append(dname)
        image_list.append(imlist)
        subject_list.append(study_list)
    return image_list, subject_list, subject_name_list

def load_timepoints(bdir, target_shape, tp_name='time??.mgz', dthresh=-1, ndilations=0):
    il, sl, sn = load_serial_cxr(bdir)

    vol_list = []
    seg_list = []
    dtrans_list = []
    for sno, ilist in enumerate(il):
        if len(ilist)>=2:
            date_list = []
            time_list = []
            for ino, im in enumerate(ilist):
                date_list.append(int(im.StudyDate))
                if hasattr(im, 'SeriesTime'):
                    time_list.append(int(im.SeriesTime.split('.')[0]))
                else:
                    time_list.append(int(im.StudyTime.split('.')[0]))

            # sort input time points by acquistion date
            ind = np.array(date_list).argsort()
            date_list2 = []
            time_list2 = []
            ilist2 = []
            for i in ind.tolist():
                date_list2.append(date_list[i])
                time_list2.append(time_list[i])
                ilist2.append(ilist[i])

            vlist = []
            slist = []
            dlist = []
            bad = False
            for ino, im in enumerate(ilist2):
                tokens = sl[sno][ind[ino]].split('/')
                fname = '/'.join(tokens[0:-2]) + '/time%2.2d.mgz' % ino
                vol = sf.load_slice(fname)
                zoomx = target_shape[0]/vol.shape[0]
                zoomy = target_shape[1]/vol.shape[1]
                vol.data = zoom(vol.data,(zoomx, zoomy),order=1)
                vlist.append(vol)

                fname = '/'.join(tokens[0:-2]) + '/time%2.2d.seg.mgz' % ino
                if os.path.exists(fname) == False:
                    print('%s missing' % fname)
                    dvol = None
                    svol = None
                    bad = True
                else:
                    svol = sf.load_slice(fname)
                    
                    svol.data = zoom(svol.data,(zoomx, zoomy),order=0)
                    u = np.unique(svol.data)
                    
                    # dilate input labels if specified by caller
                    if ndilations > 0:
                        dil_vol = np.zeros(svol.shape)
                        for l in list(u):
                            if l == 0:
                                continue
                            tmp = morph.binary_dilation(svol.data==l, iterations=ndilations)
                            dil_vol = dil_vol + l*tmp
                        svol.data = dil_vol

                    # build multiframe distance transform volume
                    dframes = []
                    for l in list(u):
                        if l == 0:
                            continue
                        dtrans = ndimage.distance_transform_edt(np.logical_not(svol.data == l))
                        if dthresh >= 0:
                            dtrans[dtrans>dthresh] = dthresh
                        dframes.append(dtrans)
                    dvol = np.transpose(np.array(dframes), (1,2,0))

                slist.append(svol)
                dlist.append(dvol)
            if bad == True:
                continue
            vol_list.append(vlist)
            seg_list.append(slist)
            dtrans_list.append(dlist)
    return vol_list, seg_list, dtrans_list, il, sl, sn


