import os
import subprocess
import glob
from nipype.pipeline.engine import Workflow, Node
from nipype.interfaces.freesurfer.preprocess import MRIConvert
from nipype.interfaces.ants import N4BiasFieldCorrection

from nipype.interfaces.ants.base import ANTSCommand, ANTSCommandInputSpec
from nipype.interfaces.base import File, traits, TraitedSpec
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import SelectFiles, DataSink

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


def preprocess(input_file, output_dir, conform=True, bias_correct=True, skullstrip=True):

    preprocess_flow = Workflow(name='preprocess', base_dir=output_dir)

    conform = Node(MRIConvert(conform=True, out_type='niigz', out_file='conformed.nii.gz'), name='conform')
    n4 = Node(N4BiasFieldCorrection(dimension=3, bspline_fitting_distance=300, shrink_factor=3,
                                    n_iterations=[50,50,30,20], output_image='n4.nii.gz' ), name='n4')
    robex = Node(ROBEX(seed=1729, stripped_image='brain.nii.gz'), name='robex')

    preprocess_flow.connect([(conform, n4, [('out_file', 'input_image')]),
                             (n4, robex, [('output_image', 'input_image')])
                             ])

    preprocess_flow.write_graph(graph2use='orig')
    conform.inputs.in_file = input_file
    preprocess_flow.run('MultiProc', plugin_args={'n_procs': 5})



# t2dir = '/autofs/space/bhim_001/users/aj660/PSACNN/data/IXI/T2/'
# t2_file_list = sorted(glob.glob(os.path.join(t2dir,'*T2.nii.gz')))
# t2_subj_list= list()
# for t2file in t2_file_list:
#     t2loc = t2file.split("/")
#     t2filename = t2loc[-1]
#     # split on T1
#     t2floc = t2filename.split("-T2")
#     t2_subj_id = t2floc[0]
#     t2_subj_list.append(t2_subj_id)
#     subprocess.call(['mkdir', '-p', os.path.join(t2dir, t2_subj_id)])
# #     # move t1_file into the subject folder
#     subprocess.call(['mv', t2file, os.path.join(t2dir, t2_subj_id)])
# #



if __name__ == "__main__":

    t2dir = '/autofs/space/bhim_001/users/aj660/PSACNN/data/IXI/T2/'
    output_dir = os.path.join(t2dir, 'processed')
    subprocess.call(['mkdir', '-p', output_dir])
    t2_file_list = sorted(glob.glob(os.path.join(t2dir, 'IXI*')))


    # t1_file_list = sorted(glob.glob(os.path.join(t1dir,'*T1.nii.gz')))
    t2_subj_list= list()
    for t2file in t2_file_list:
        t2loc = t2file.split("/")
        t2_subj_id = t2loc[-1]
        # split on T1
    #     t1floc = t1filename.split("-T1")
    #     t1_subj_id = t1floc[0]
        t2_subj_list.append(t2_subj_id)
    #     subprocess.call(['mkdir', '-p', os.path.join(t1dir, t1_subj_id)])
    #     # move t1_file into the subject folder
    #     subprocess.call(['mv', t1file, os.path.join(t1dir, t1_subj_id)])
    #


    # Infosource - a function free node to iterate over the list of subject names
    infosource = Node(IdentityInterface(fields=['subject_id']),
                      name="infosource")
    infosource.iterables = [('subject_id', t2_subj_list)]


    # SelectFiles - to grab the data (alternativ to DataGrabber)
    anat_file = os.path.join('{subject_id}', '{subject_id}-T2.nii.gz')
    templates = {'anat': anat_file}

    selectfiles = Node(SelectFiles(templates,
                                   base_directory=t2dir),
                       name="selectfiles")

    datasink = Node(DataSink(base_directory=t2dir,
                             container=output_dir),
                    name="datasink")
    substitutions = [('_subject_id_', '')]
    datasink.inputs.substitutions = substitutions


    # t2dir = '/autofs/space/bhim_001/users/aj660/PSACNN/data/IXI/T2/'
    # t2_file_list = sorted(glob.glob(os.path.join(t1dir, '*T2.nii.gz')))
    #
    # input_file = '/autofs/space/bhim_001/users/aj660/PSACNN/data/IXI/T2/IXI002-Guys-0828-T2.nii.gz'
    # output_dir = '/autofs/space/bhim_001/users/aj660/psacnn_brain_segmentation/test_output/IXI002-Guys-0828-T2'
    # subprocess.call(['mkdir', '-p', output_dir])

    preprocess_flow = Workflow(name='preprocess', base_dir=output_dir)

    conform = Node(MRIConvert(conform=True, out_type='niigz', out_file='conformed.nii.gz'), name='conform')
    n4 = Node(N4BiasFieldCorrection(dimension=3, bspline_fitting_distance=300, shrink_factor=3,
                                    n_iterations=[50,50,30,20], output_image='n4.nii.gz' ), name='n4')
    robex = Node(ROBEX(seed=1729, stripped_image='brain.nii.gz'), name='robex')

    preprocess_flow.connect([(infosource, selectfiles, [('subject_id', 'subject_id')]),
                             (selectfiles, conform, [('anat', 'in_file')]),
                             (conform, n4, [('out_file', 'input_image')]),
                             (n4, robex, [('output_image', 'input_image')]),
                             (robex, datasink, [('stripped_image', 'preprocess.brain'),
                                                ('mask_image', 'preprocess.mask')])
                             ])

    preprocess_flow.write_graph(graph2use='colored', format='png', simple_form=True)
    # conform.inputs.in_file = input_file
    preprocess_flow.run('MultiProc', plugin_args={'n_procs': 20})



