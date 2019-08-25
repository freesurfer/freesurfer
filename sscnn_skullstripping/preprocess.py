import numpy as np
import subprocess
from os.path import join as opj

from sscnn_skullstripping.deeplearn_utils import DeepImageSynth

import os
import subprocess
from nipype.pipeline.engine import Workflow, Node
from nipype.interfaces.freesurfer.preprocess import MRIConvert

from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec
from nipype.interfaces.base import File, traits, TraitedSpec, Directory
import nibabel as nib

import collections
from skimage.measure import label

def prob_to_hard(cor_img_data, sag_img_data, ax_img_data, p_cor_v, p_sag_v, p_ax_v):
    cor_ax_img_data = np.transpose(ax_img_data[:, :, :, 1], (0, 2, 1))
    cor_img_data = cor_img_data[:, :, :, 1]
    cor_sag_img_data = np.transpose(sag_img_data[:, :, :, 1], (2, 1, 0))


    combine_img_data = p_cor_v * cor_img_data + \
                       p_sag_v * cor_sag_img_data + \
                       p_ax_v * cor_ax_img_data

    hard_img_data = np.zeros(combine_img_data.shape)
    hard_img_data[combine_img_data > 0.5] = 1
    return hard_img_data


def predict_segmentation(input_file, output_dir, contrast = 't1w',
                         ax_model_file=None, cor_model_file=None, sag_model_file=None,
                               save_label_image=True, save_prob_image=False,
                        batch_size=4):

    brain_labels = [0,1] #np.loadtxt('aseg_labels.txt')
    num_input_channels = 1
    channel_names = [contrast]
    feature_shape = (256, 256, num_input_channels)

    filter_size = (7, 7)
    f = 32
    d = 6




    model_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'model_files')
    if ax_model_file is None:
        ax_model_file = opj(model_dir, 'ax_sscnn.h5')
    if cor_model_file is None:
        cor_model_file = opj(model_dir, 'cor_sscnn.h5')
    if sag_model_file is None:
        sag_model_file = opj(model_dir, 'sag_sscnn.h5')

    loss = 'dice_coef_loss2'
    ax_curr_unet = DeepImageSynth.DeepImageSynth.from_file(ax_model_file,loss, net='unet_2d_v1',
                                                           n_labels=len(brain_labels),labels=brain_labels,
                                                           storage_loc='disk',  temp_folder=model_dir,
                                                           rob_standardize=True,
                                                           wmp_standardize=False)

    cor_curr_unet = DeepImageSynth.DeepImageSynth.from_file(cor_model_file,loss, net='unet_2d_v1',
                                                            n_labels=len(brain_labels),labels=brain_labels,
                                                            storage_loc='disk',  temp_folder=model_dir,
                                                            rob_standardize=True,
                                                            wmp_standardize=False)

    sag_curr_unet = DeepImageSynth.DeepImageSynth.from_file(sag_model_file,loss, net='unet_2d_v1',
                                                            n_labels=len(brain_labels),labels=brain_labels,
                                                            storage_loc='disk',  temp_folder=model_dir,
                                                            rob_standardize=True,
                                                            wmp_standardize=False)



    ax_curr_unet.rob_standardize = True
    ax_curr_unet.wmp_standardize = False

    cor_curr_unet.rob_standardize = True
    cor_curr_unet.wmp_standardize = False

    sag_curr_unet.rob_standardize = True
    sag_curr_unet.wmp_standardize = False

    ax_out_membership_file = opj(output_dir, "sscnn_ax_prob.mgz")
    ax_out_hard_file = opj(output_dir, "sscnn_ax_label.mgz")
    cor_out_membership_file = opj(output_dir, "sscnn_cor_prob.mgz")
    cor_out_hard_file = opj(output_dir, "sscnn_cor_label.mgz")
    sag_out_membership_file = opj(output_dir, "sscnn_sag_prob.mgz")
    sag_out_hard_file = opj(output_dir, "sscnn_sag_label.mgz")

    # loss = 'dice_coef_loss2'
    # curr_unet = DeepImageSynth.DeepImageSynth.from_file(model_file, loss, n_labels=len(aseg_labels), labels=aseg_labels,
    #                                                     storage_loc='disk', temp_folder=model_dir,
    #                                                     num_input_channels=num_channels, channel_names=channel_names,
    #                                                     out_channel_names=['seg'], initial_learning_rate=0.00001,
    #                                                     use_patches=True,
    #                                                     wmp_standardize=True, rob_standardize=True, fcn=True,
    #                                                     num_gpus=1,
    #                                                     preprocessing=False, augment=False, num_outputs=1,
    #                                                     nmr_augment=True, use_tal=False, rare_label_list=[]
    #

    #                                                      )

    test_files = list()
    for iter_channel in range(num_input_channels):
        test_files.append(input_file)

    ax_curr_unet.predict_slice_segmentation(test_files, channel_names, 'axial', ax_out_membership_file,
                                            ax_out_hard_file)
    cor_curr_unet.predict_slice_segmentation(test_files, channel_names, 'coronal', cor_out_membership_file,
                                             cor_out_hard_file)
    sag_curr_unet.predict_slice_segmentation(test_files, channel_names, 'sagittal', sag_out_membership_file,
                                             sag_out_hard_file)

    # read in hard seg files

    p_c = 0.44
    p_s = 0.33
    p_a = 0.23

    ax_img = nib.load(ax_out_membership_file)
    ax_img_data = ax_img.get_data()

    cor_img = nib.load(cor_out_membership_file)
    cor_img_data = cor_img.get_data()

    sag_img = nib.load(sag_out_membership_file)
    sag_img_data = sag_img.get_data()

    hard_img_data = prob_to_hard(cor_img_data, sag_img_data, ax_img_data, p_c, p_s, p_a)

    # add a post processing function to only choose the largest connected component
    label_img_data = label(np.int8(hard_img_data), neighbors=8)
    freq = collections.Counter(label_img_data.flatten())


    l_idx = np.argmax(list(freq.values())[1:])
    big_label = list(freq.keys())[1:][l_idx]
    hard_img_data[label_img_data != big_label] = 0

    # out_img = nib.MGHImage(hard_img_data, cor_img.affine, cor_img.header)

    # majority voting. can be better.

    print('creating output directory: ' + output_dir)
    subprocess.call(['mkdir', '-p', output_dir])

    # out_hard_file = opj(output_dir, "sscnn_strip.mgz")
    # nib.save(out_img, out_hard_file)
    return hard_img_data


class SSCNNInputSpec(BaseInterfaceInputSpec):
    import os
    input_image = File(exists=True, desc='intensity image to be segmented', argstr='-i %s', position=0, mandatory=True)
    output_dir = Directory(exists=True, desc='output directory to store results', argstr='-o %s', position=1, mandatory=True)
    contrast = traits.String(desc='contrast of input image. t1w or t2w', argstr='-c %s', position=2, mandatory=True)

    save_label_image = traits.Bool(True, desc='saves label image', argstr='-ol', position=5,
                                   mandatory=False)

    save_prob_image = traits.Bool(False, desc='saves probability map. can be 1.5GB', argstr='-os', position=6,
                                   mandatory=False)

    # model_file = File(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'model_files', 'flt1t2_psacnn_p96.h5'),
    #                   exists=True, desc='model file', argstr='-m %s', position=7, mandatory=False, usedefault=True)

    batch_size = traits.Int(4, desc='batch size for GPU processing', argstr='-batch_size %d', position=8, usedefault=True)


    # label_list_file = File(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'model_files','aseg_labels_skullstripped.txt')
    #                        , exists=True, desc='label list file', argstr='-l %s', position=9, mandatory=False, usedefault=True)
    #
    # sample_rate = traits.Int(20000, desc='predict only for N/sample_rate number of patches', argstr='-sample_rate %d', position=10, usedefault=True)


class SSCNNOutputSpec(TraitedSpec):
    label_image = File(exists=True, desc='hard segmentation label image')
    # prob_image = File(exists=True, desc='probability map image')


class SSCNN(BaseInterface):

    """
    Runs SSCNN via the sscnn_segment function.

    Examples:
    ---------
    'run psacnn_segment.py -i t1.nii.gz -o output_dir/psacnn -gpu 0 -c 't2w' -p 96
   '
    >>> sscnn = SSCNN()
    >>> sscnn.inputs.input_image = 't1.nii.gz'
    >>> sscnn.inputs. = 'output_dir'
    >>> sscnn.inputs.contrast = 't1w'
    >>> sscnn.cmdline
    """

    input_spec = SSCNNInputSpec
    output_spec = SSCNNOutputSpec

    def _run_interface(self, runtime):


        '''(input_file, output_dir, contrast = 't1w',
                         ax_model_file=None, cor_model_file=None, sag_model_file=None,
                               save_label_image=True, save_prob_image=False,
                        batch_size=4)'''
        label_image = predict_segmentation(self.inputs.input_image,
                                           self.inputs.output_dir,
                                           contrast=self.inputs.contrast,
                                           ax_model_file=None,
                                           cor_model_file=None,
                                           sag_model_file=None,
                                           save_label_image=True,
                                           save_prob_image=True,
                                           batch_size=4)


        in_img = nib.load(self.inputs.input_image)
        out_img = nib.Nifti1Image(label_image, in_img.affine, in_img.header)
        out_file = os.path.join(self.inputs.output_dir, 'sscnn_skullstrip.nii.gz')
        nib.save(out_img, out_file)




        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        base_dir = self.inputs.output_dir
        outputs["label_image"] = os.path.join(base_dir, 'sscnn_skullstrip.nii.gz')
        return outputs



def sscnn_workflow(input_file, output_dir, contrast='t1w', use_gpu=True,
                   gpu_id=0, save_label_image=False, save_prob_image=False, batch_size=4):

    subprocess.call(['mkdir', '-p', output_dir])
    if use_gpu == False:
        os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        os.environ["CUDA_VISIBLE_DEVICES"] = "-1"
        gpu_id = -1
        batch_size = 16
    else:
        os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
        batch_size = 4


    preprocess_flow = Workflow(name='preprocess', base_dir=output_dir)
    conform = Node(MRIConvert(conform=True, out_type='niigz', out_file='conformed.nii.gz'), name='conform')
    sscnn = Node(SSCNN(output_dir=output_dir, contrast=contrast, batch_size=batch_size,
                         save_label_image=save_label_image, save_prob_image=save_prob_image), name='sscnn')

    preprocess_flow.connect([(conform, sscnn, [('out_file', 'input_image')]),])

    preprocess_flow.write_graph(graph2use='orig')
    conform.inputs.in_file = input_file
    preprocess_flow.run('MultiProc', plugin_args={'n_procs': 16})






if __name__ == "__main__":



    # input_file = '/autofs/space/bhim_001/users/aj660/PSACNN/data/IXI/T2/preprocess/brain/IXI511-HH-2238/brain.nii.gz'
    # input_file = '/autofs/space/bhim_001/users/aj660/PSACNN/data/IXI/T2/processed/preprocess/brain/IXI511-HH-2238/brain.nii.gz'
    input_file = '/autofs/space/vault_007/users/lzollei/AmodJog/data/NadineGaab-Bangladesh/01/mprage.nii.gz'
    output_dir = '/autofs/space/bhim_001/users/aj660/tmp'
    subprocess.call(['mkdir', '-p', output_dir])

    sscnn_workflow(input_file, output_dir,  contrast='t1w', use_gpu=False,
                    gpu_id=0, save_label_image=False, save_prob_image=False,batch_size=4)
