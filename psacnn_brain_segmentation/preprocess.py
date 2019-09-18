import numpy as np
import subprocess
from os.path import join as opj

from psacnn_brain_segmentation.deeplearn_utils import DeepImageSynth

import os
import subprocess
from nipype.pipeline.engine import Workflow, Node
from nipype.interfaces.freesurfer.preprocess import MRIConvert
from nipype.interfaces.ants import N4BiasFieldCorrection

from nipype.interfaces.ants.base import ANTSCommand, ANTSCommandInputSpec
from nipype.interfaces.base import BaseInterface, BaseInterfaceInputSpec
from nipype.interfaces.base import File, traits, TraitedSpec, Directory
import nibabel as nib



def predict_segmentation(input_file, output_dir, contrast = 't1w',
                         model_file=None, save_label_image=True, save_prob_image=False,
                         patch_dim=96, batch_size=4, label_list_file=None, sample_rate=20000):

    aseg_labels = np.loadtxt(label_list_file)
    num_channels = 1
    channel_names = [contrast]
    feature_shape = (96, 96, 96, num_channels)
    f = 32
    d = 5



    model_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'model_files')
    if model_file is None:
       model_file = opj(model_dir, 'flt1t2_psacnn_p'+str(patch_dim)+'.h5')

    loss = 'dice_coef_loss2'
    curr_unet = DeepImageSynth.DeepImageSynth.from_file(model_file, loss, net='unet_2d_v1',n_labels=len(aseg_labels), labels=aseg_labels,
                                                        storage_loc='disk', temp_folder=model_dir,
                                                        num_input_channels=num_channels, channel_names=channel_names,
                                                        out_channel_names=['seg'], initial_learning_rate=0.00001,
                                                        use_patches=True,
                                                        wmp_standardize=True, rob_standardize=True, fcn=True,
                                                        num_gpus=1,
                                                        preprocessing=False, augment=False, num_outputs=1,
                                                        nmr_augment=False, use_tal=False, rare_label_list=[]
                                                        )

    print('creating output directory: ' + output_dir)
    subprocess.call(['mkdir', '-p', output_dir])

    test_files = list()
    if num_channels == 1:
        test_files.append(input_file)

    out_membership_file = opj(output_dir, "psacnn_soft.mgz")
    out_hard_file = opj(output_dir, "psacnnseg.mgz")
    test_channel_names = list()
    test_channel_names.append(contrast)
    label_image = curr_unet.predict_segmentation(test_files, test_channel_names, out_membership_file, out_hard_file,
                                   step_size=[16, 16, 16], batch_size=batch_size, center_voxel=True, sampling_rate=sample_rate,
                                   save_label_image=save_label_image, save_prob_image=save_prob_image)


    return label_image

# def run_segmentation_pipeline(input_file, output_dir, use_preprocess=False,
#                               contrast='t1w', gpu=0, patch_size=96,model_file=None, output_soft=False):


# Skull Stripping
class ROBEXInputSpec(ANTSCommandInputSpec):
    input_image = File(exists=True, desc='T1 image to be skull stripped', argstr='%s', position=0, mandatory=True)
    stripped_image = File(desc='output file name of stripped image', name_source=['input_image'], hash_files=False,
                          keep_extension=True, name_template='%s_stripped', argstr='%s', position=1)
    mask_image = File(desc='output file name of binary mask', name_source=['input_image'], hash_files=False,
                      keep_extension=True, name_template='%s_mask', argstr='%s', position=2)
    seed = traits.Int(desc='seed value for non-deterministic functions', argstr='%d', position=3)


class ROBEXOutputSpec(TraitedSpec):
    stripped_image = File(exists=True, desc='skull stripped image')
    mask_image = File(exists=True, desc='binary brain mask')


class ROBEX(ANTSCommand):
    """
    Runs ROBEX via the runROBEX.sh script included in the ROBEX package.

    Examples:
    ---------



    >>> robex = ROBEX()
    >>> robex.inputs.input_image = 't1.nii.gz'
    >>> robex.cmdline
    'runROBEX.sh t1.nii.gz t1_stripped.nii.gz t1_mask.nii.gz'
    """

    input_spec = ROBEXInputSpec
    output_spec = ROBEXOutputSpec
    _cmd = 'runROBEX.sh'



class PSACNNInputSpec(BaseInterfaceInputSpec):
    import os
    input_image = File(exists=True, desc='intensity image to be segmented', argstr='-i %s', position=0, mandatory=True)
    output_dir = Directory(exists=True, desc='output directory to store results', argstr='-o %s', position=1, mandatory=True)
    contrast = traits.String(desc='contrast of input image. t1w or t2w', argstr='-c %s', position=2, mandatory=True)
    # gpu = traits.Int(desc='seed value for non-deterministic functions', argstr='-gpu %d', position=3)
    patch_size = traits.Enum(96, 64, desc='dimension of cube patch one of 64 or 96. Default 96', argstr='-p %d', position=4, usedefault=True)

    save_label_image = traits.Bool(True, desc='saves label image', argstr='-ol', position=5,
                                   mandatory=False)

    save_prob_image = traits.Bool(False, desc='saves probability map. can be 1.5GB', argstr='-os', position=6,
                                   mandatory=False)

    model_file = File(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'model_files', 't1t2_psacnn_p96_v1.h5'),
                      exists=True, desc='model file', argstr='-m %s', position=7, mandatory=False, usedefault=True)

    batch_size = traits.Int(4, desc='batch size for GPU processing', argstr='-batch_size %d', position=8, usedefault=True)


    label_list_file = File(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'model_files','aseg_labels_skullstripped.txt')
                           , exists=True, desc='label list file', argstr='-l %s', position=9, mandatory=False, usedefault=True)

    sample_rate = traits.Int(20000, desc='predict only for N/sample_rate number of patches', argstr='-sample_rate %d', position=10, usedefault=True)


class PSACNNOutputSpec(TraitedSpec):
    label_image = File(exists=True, desc='hard segmentation label image')
    # prob_image = File(exists=True, desc='probability map image')


class PSACNN(BaseInterface):

    """
    Runs PSACNN via the psacnn_segment function.

    Examples:
    ---------
    'run psacnn_segment.py -i t1.nii.gz -o output_dir/psacnn -gpu 0 -c 't2w' -p 96
   '
    >>> psacnn = PSACNN()
    >>> psacnn.inputs.input_image = 't1.nii.gz'
    >>> psacnn.inputs. = 'output_dir'
    >>> psacnn.inputs.contrast = 't1w'
    >>> robex.cmdline
    'runROBEX.sh t1.nii.gz t1_stripped.nii.gz t1_mask.nii.gz'
    """

    input_spec = PSACNNInputSpec
    output_spec = PSACNNOutputSpec

    def _run_interface(self, runtime):



        label_image = predict_segmentation(self.inputs.input_image, self.inputs.output_dir, contrast=self.inputs.contrast,
                             model_file=self.inputs.model_file, save_label_image=self.inputs.save_label_image,
                             save_prob_image=self.inputs.save_prob_image,patch_dim=self.inputs.patch_size,
                             batch_size=self.inputs.batch_size, label_list_file=self.inputs.label_list_file,
                                           sample_rate=self.inputs.sample_rate)


        in_img = nib.load(self.inputs.input_image)
        out_img = nib.Nifti1Image(label_image, in_img.affine, in_img.header)
        out_file = os.path.join(self.inputs.output_dir, 'psacnnseg.nii.gz')
        nib.save(out_img, out_file)




        return runtime

    def _list_outputs(self):
        outputs = self._outputs().get()
        base_dir = self.inputs.output_dir
        outputs["label_image"] = os.path.join(base_dir, 'psacnnseg.nii.gz')
        return outputs



def psacnn_workflow(input_file, output_dir, use_preprocess=True, model_file=None, contrast='t1w', use_gpu=True,
                    gpu_id=0, save_label_image=False, save_prob_image=False, patch_size=96, batch_size=4, sample_rate=20000):

    subprocess.call(['mkdir', '-p', output_dir])
    if use_gpu == False:
        os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        gpu_id = -1
        batch_size = 16
        sample_rate = 40000
    else:
        os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"
        os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu_id)
        batch_size = 4
        sample_rate = 20000


    if use_preprocess == True:
        preprocess_flow = Workflow(name='preprocess', base_dir=output_dir)

        conform = Node(MRIConvert(conform=True, out_type='niigz', out_file='conformed.nii.gz'), name='conform')
        n4 = Node(N4BiasFieldCorrection(dimension=3, bspline_fitting_distance=300, shrink_factor=3,
                                        n_iterations=[50,50,30,20], output_image='n4.nii.gz' ), name='n4')
        robex = Node(ROBEX(seed=1729, stripped_image='brain.nii.gz'), name='robex')

        psacnn = Node(PSACNN(output_dir=output_dir, contrast=contrast, patch_size=patch_size, batch_size=batch_size,
                             save_label_image=save_label_image, save_prob_image=save_prob_image, sample_rate=sample_rate), name='psacnn')

        preprocess_flow.connect([(conform, n4, [('out_file', 'input_image')]),
                                 (n4, robex, [('output_image', 'input_image')]),
                                 (robex, psacnn, [('stripped_image', 'input_image')])
                                 ])

        preprocess_flow.write_graph(graph2use='orig')
        conform.inputs.in_file = input_file
        preprocess_flow.run('MultiProc', plugin_args={'n_procs': 16})
    else:

        psacnn = PSACNN(input_image=input_file, output_dir=output_dir, contrast=contrast, patch_size=patch_size,
                        batch_size=batch_size, save_label_image=save_label_image,save_prob_image=save_prob_image,
                        sample_rate=sample_rate)
        # psacnn.inputs.input_image = input_file
        # psacnn.inputs.output_dir = output_dir
        # psacnn.inputs.contrast = contrast
        # psacnn.inputs.patch_size = patch_size
        # psacnn.inputs.batch_size = batch_size
        # psacnn.inputs.save_label_image = save_label_image
        # psacnn.inputs.save_prob_image = save_prob_image
        # psacnn.inputs.sample_rate = sample_rate

        psacnn.run()






if __name__ == "__main__":



    input_file = '/autofs/space/bhim_001/users/aj660/PSACNN/data/IXI/T1/preprocess/brain/IXI511-HH-2238/brain.nii.gz'
    output_dir = '/autofs/space/bhim_001/users/aj660/psacnn_brain_segmentation/test_output/preproc_T1_IXI511-HH-2238'
    subprocess.call(['mkdir', '-p', output_dir])

    psacnn_workflow(input_file, output_dir, use_preprocess=False, model_file=None, contrast='t1w', use_gpu=True,
                    gpu_id=0, save_label_image=False, save_prob_image=False, patch_size=96, batch_size=4,
                    sample_rate=20000)
