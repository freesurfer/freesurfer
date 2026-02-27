import argparse
import os
import torch
import numpy as np
import glob
import random
from photo_reconstruction.image_utils import MRIread

# Parse arguments of command line
def get_arguments():
    # Parse arguments (tons of them!)
    parser = argparse.ArgumentParser(
        description="Code for 3D photo reconstruction (Tregidgo, et al., MICCAI 2020)"
    )

    parser.add_argument("--input_photo_dir", type=str, required=True,
                        help="Directory with input photos (required)")

    parser.add_argument("--input_segmentation_dir", type=str, required=True,
                        help="Directory with input slab masks / segmentations (required)")

    parser.add_argument("--ref_mri", type=str, default=None,
                        help="Reference MRI scan (if available)")

    parser.add_argument("--ref_mri_synthseg", type=str, default=None,
                        help="SynthSeg of reference MRI scan (will be computed if it does not exist)")

    parser.add_argument("--ref_mri_synthsr", type=str, default=None,
                        help="SynthSR of reference MRI; only needed if MRI is not 1mm T1 (will be computed if it does not exist)")

    parser.add_argument("--low_field_synthsr", dest="low_field_synthsr", action="store_true",
                        help="Use low-field version of SynthSR (eg for Hyperfine scans)")
    parser.set_defaults(low_field_synthsr=False)

    parser.add_argument("--input_roi_dir", type=str, default=None,
                        help="Directory with ROi masks to deform (optional)")

    parser.add_argument("--ref_mesh", type=str, default=None,
                        help="Reference surface mesh (if available)")

    parser.add_argument("--mesh_reorient_with_indices", type=str, default=None,
                        help="Indices to reorient mesh (see PhotoTools page in FreeSurfer wiki)")

    parser.add_argument("--fresh_tissue", dest="fresh_tissue", action="store_true",
                        help="Uses more lenient regularizers to accommodate fresh tissue")
    parser.set_defaults(fresth_tissue=False)

    parser.add_argument("--photos_of_posterior_side", dest="posterior_side", action="store_true",
                        help="Use when photos are taken of posterior side of slabs (default is anterior side)")
    parser.set_defaults(posterior_side=False)

    parser.add_argument("--order_posterior_to_anterior", dest="posterior_to_anterior", action="store_true",
                        help="Use when photos are ordered from posterior to anterior (default is anterior to posterior)")
    parser.set_defaults(posterior_to_anterior=False)

    parser.add_argument("--hemisphere", type=str, required=True,
                        help="hemisphere; must be left, right, or both (required)")

    parser.add_argument("--slice_thickness", type=float, required=True,
                        help="Slice thickness in mm (required); will be finetuned if possible")

    parser.add_argument("--photo_resolution", type=float, required=True,
                        help="Resolution of the photos in mm (required)")

    parser.add_argument("--initial_stretch_factor_lr_photos", type=float, default=1.0,
                        help="Initialize stretch of photos in left-right direction by this factor.")

    parser.add_argument("--no_z_stretch", dest="no_z_stretch", action="store_true",
                        help="Use when you are certain of slice thickness and/or photos are outside the mesh")
    parser.set_defaults(no_z_stretch=False)

    parser.add_argument("--use_svf_for_photos", dest="use_svf_for_photos", action="store_true",
                        help="Uses a diffeomorphic deformation model (SVF) for the photos")
    parser.set_defaults(use_svf_for_photos=False)

    parser.add_argument("--stretch_factor_lr_mesh", type=float, default=1.0,
                        help="Stretch mesh in left-right direction by this factor.")

    parser.add_argument("--weights", type=str, default=None,
                        help="CSV file with slab weights")

    parser.add_argument("--thickness_cap", type=float, default=0,
                        help="Maximum thickness that you allow when estimating thicknesses from weights")

    parser.add_argument("--output_directory", type=str, required=True,
                        help="Output directory with reconstructed photo volume and reference (required)")

    parser.add_argument("--equalize_images", dest="equalize_images", action="store_true",
                        help="Use to equalize images (useful if they have low contrast)",
                        )
    parser.set_defaults(equalize_images=False)

    parser.add_argument("--skip_bfgs", dest="skip_bfgs", action="store_true",
                        help="Use to skip BFGS finetuning (after optimization with Adam)",
                        )
    parser.set_defaults(skip_bfgs=False)

    parser.add_argument("--threads", type=int, default=-1,
                        help="Number of cores to be used. Default is 1. You can use -1 to use all available cores")

    parser.add_argument("--gpu", type=int, default=None,
                        help="Index of GPU to use (default is None)")

    # This is the option to deform the surfaces / segmentations from the reference's FreeSurfer directory
    # (which obviously required running the reference through FreeSurfer first!).
    parser.add_argument("--deform_recon_dir", type=str, default=None,
                        help="Directory with FS dir of reference, to deform surfaces etc (expects to find deform_recon_dir/surf) (optional)")

    # Control point spacing for nonlinear deformation models
    parser.add_argument("--cp_spacing_2d", type=float, default=None, help="(Advanced) Control point spacing for 2D deformation")
    parser.add_argument("--cp_spacing_3d", type=float, default=None, help="(Advanced) Control point spacing for 3D deformation")

    # weights of different terms in loss
    parser.add_argument("--k_lncc_mri", type=float, default=None,
                        help="(Advanced) Weight of LNCC between reference MRI and reconstruction")
    parser.add_argument("--k_dice_mri", type=float, default=None,
                        help="(Advanced) Weight of Dice between masks of reference and reconstruction")
    parser.add_argument("--k_dif_slice_loss", type=float, default=None,
                        help="(Advanced) Weight of SSD between consecutive slices of reconstruction")
    parser.add_argument("--k_mesh_loss", type=float, default=None,
                        help="(Advanced) Weight of absolute distance to edge of masks from mesh vertices")
    parser.add_argument("--k_regularizer", type=float, default=None,
                        help="(Advanced) Weight of regularizer of log_det(affine matrices)")
    parser.add_argument("--k_regularizer_nonlin", type=float, default=None,
                        help="(Advanced) Weight of regularizer of 2D nonlinear deformation of photos")
    parser.add_argument("--k_regularizer_nonlin3d", type=float, default=None,
                        help="(Advanced) Weight of regularizer of 3D nonlinear deformatin of reference")
    parser.add_argument("--k_regularizer_sz", type=float, default=None,
                        help="(Advanced) Weight of regularizer of stretch in AP direction")

    arguments = parser.parse_args()

    # Some checks ...
    if (arguments.ref_mri_synthseg is not None) and (arguments.ref_mri is None):
        raise Exception('You provided an MRI segmentation but no MRI scan')

    if (arguments.ref_mri_synthsr is not None) and (arguments.ref_mri is None) and (not os.path.exists(ref_mri)):
        raise Exception('You provided an MRI synthsr that does not exist but no MRI scan to compute it from')

    if (arguments.mesh_reorient_with_indices is not None) and (arguments.ref_mesh is None):
        raise Exception('You provided mesh indices but no mesh')

    if (arguments.mesh_reorient_with_indices is None) and (arguments.ref_mesh is not None):
        raise Exception('You provided mesh but no indices to reorient it (currently not supported)')

    if (arguments.hemisphere!='left') and (arguments.hemisphere!='right') and (arguments.hemisphere!='both'):
        raise Exception('Hemisphere must be left, right, or both')

    if (arguments.ref_mri is not None) and (not os.path.exists(arguments.ref_mri)):
        raise Exception('Input MRI does not exist')

    return arguments

# Configure devices and number of CPU threads
def configure_gpu_and_cpu(gpu, threads):
    if gpu is None:
        print("Using the CPU")
        os.environ["CUDA_VISIBLE_DEVICES"] = ""
        device = torch.device("cpu")
    else:
        os.environ["CUDA_VISIBLE_DEVICES"] = str(gpu)
        if torch.cuda.is_available():
            print("Using GPU device " + str(gpu))
            device = torch.device("cuda:0")
        else:
            print(
                "Tried to use GPU device "
                + str(gpu)
                + " but failed; using CPU instead"
            )
            device = torch.device("cpu")

    if threads == 1:
        print('using 1 thread')
    elif threads<0:
        threads = os.cpu_count()
        print('using all available threads ( %s )' % threads)
    else:
        print('using %s threads' % threads)
    torch.set_num_threads(threads)

    return device

# Dumb function to figure out the prefix we need to source freesurfer in system calls
def get_fs_prefix():
    print('Quickly detecting shell; trying bash')
    fshome = os.getenv('FREESURFER_HOME')
    fsprefix = 'export FREESURFER_HOME=' + fshome + '; source ' + fshome + '/SetUpFreeSurfer.sh ; '
    a = os.system(fsprefix + ' mri_convert --help  >/dev/null')
    if a > 0:
        print('bash failed; trying with tcsh...')
        fsprefix = 'setenv FREESURFER_HOME ' + fshome + '; source ' + fshome + '/SetUpFreeSurfer.csh ; '
        a = os.system(fsprefix + ' mri_convert --help  >/dev/null')
        if a > 0:
            raise Exception('both shells failed; exitting')
    return fsprefix

# Set all seeds, for reproducibility purposes
def seed_all(seed):
    # https://discuss.pytorch.org/t/reproducibility-with-all-the-bells-and-whistles/81097
    seed = 0 if not seed else seed
    print("[ Using Seed : ", seed, " ]")

    torch.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.cuda.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


# Get list of files with photos and masks, in the correct order
def get_photo_and_seg_lists(input_photo_dir, input_segmentation_dir):

    d_i = glob.glob(input_photo_dir + "/*.jpg")
    if len(d_i) == 0:
        d_i = glob.glob(input_photo_dir + "/*.jpeg")
    if len(d_i) == 0:
        d_i = glob.glob(input_photo_dir + "/*.tif")
    if len(d_i) == 0:
        d_i = glob.glob(input_photo_dir + "/*.tiff")
    if len(d_i) == 0:
        d_i = glob.glob(input_photo_dir + "/*.png")
    if len(d_i) == 0:
        d_i = glob.glob(input_photo_dir + "/*.JPG")
    if len(d_i) == 0:
        d_i = glob.glob(input_photo_dir + "/*.JPEG")
    if len(d_i) == 0:
        d_i = glob.glob(input_photo_dir + "/*.TIF")
    if len(d_i) == 0:
        d_i = glob.glob(input_photo_dir + "/*.TIFF")
    if len(d_i) == 0:
        d_i = glob.glob(input_photo_dir + "/*.PNG")
    d_i = sorted(d_i)

    d_s = glob.glob(input_segmentation_dir + "/*.mat")  # try a bunch of extensions
    if len(d_s) == 0:
        d_s = glob.glob(input_segmentation_dir + "/*.npy")
    if len(d_s) == 0:
        d_s = glob.glob(input_segmentation_dir + "/*.MAT")
    if len(d_s) == 0:
        d_s = glob.glob(input_segmentation_dir + "/*.NPY")
    d_s = sorted(d_s)

    return d_i, d_s

# Prepares reference volumes, runs SynthSR / SynthSeg if needed, etc
def prepare_reference_volumes(ref_mri, ref_mri_synthsr, low_field_synthsr, ref_mri_synthseg, output_directory, fsprefix, threads):
    print('Preparing 3D MRI reference ...')

    # Run SynthSR if needed
    if ref_mri_synthsr is None:
        input_scan = ref_mri
        print('  ***IMPORTANT***   No SynthSR provided; reference had better be a 1mm T1 scan (e.g., MP-RAGE) or pre-computed SynthSR of MNI template')
    else:
        input_scan = ref_mri_synthsr
        if os.path.exists(input_scan):
            print('  SynthSR already exists; no need to compute it')
        else:
            print('  Running SynhtSR on input scan')
            cmd = fsprefix + ' mri_synthsr --i ' + ref_mri + ' --o ' + input_scan + ' --threads ' + str(threads)
            if low_field_synthsr:
                cmd += ' --lowfield '
            a = os.system(cmd + ' >/dev/null')
            if a > 0:
                raise Exception('mri_synthsr failed; exiting')

    # Run SynthSeg if needed
    if ref_mri_synthseg is None:
        input_seg = output_directory + '/synthseg.nii.gz'
        print('  SynthSeg not provided; defaulting to: ' + input_seg)
    else:
        input_seg = ref_mri_synthseg

    resampled_seg = output_directory + '/synthseg.resampled.nii.gz'
    if os.path.exists(input_seg):
        print('  SynthSeg already exists; no need to compute it')
    else:
        print('  Running SynhtSeg on reference volume')
        cmd = fsprefix + ' mri_synthseg --i ' + input_scan + ' --o ' + input_seg + ' --robust --threads ' + str(threads)
        a = os.system(cmd + ' >/dev/null')
        if a > 0:
            raise Exception('mri_synthseg failed; exiting')
    # Resample to space of reference volume
    cmd = fsprefix + ' mri_convert ' + input_seg + ' ' + resampled_seg + ' -rl ' + input_scan + ' -rt nearest -odt float'
    a = os.system(cmd + ' >/dev/null')
    if a > 0:
        raise Exception('mri_convert failed; exiting')

    REF, REFaff = MRIread(input_scan)
    if np.isnan(REF).any():
        print("There are NaNs is the reference volume; we replaced them by zeros")
        REF[np.isnan(REF)] = 0
    REF = np.squeeze(REF)
    REF_SSEG, _ = MRIread(resampled_seg)
    os.system('rm -rf ' + resampled_seg + ' >/dev/null')
    
    return REF, REF_SSEG, REFaff

# Comes up with values for regularizers, contros points spacings, etc, depending on available references, fresh vs fixed tissue, etc
def adjust_settings(arguments):

    print('Adjusting settings automatically')

    mri_present = (arguments.ref_mri is not None) or (arguments.ref_mri_synthsr is not None)
    mesh_present = (arguments.ref_mesh is not None)
    fresh = arguments.fresh_tissue

    # Loss weights are constant
    k_lncc_mri = k_dice_mri = k_mesh_loss = 1.0
    k_dif_slice_loss = 3.0

    # Affine regularizer: 0.1 seems pretty good
    k_regularizer = 0.1

    # A-P stretch regularizer: after many tests, 0.0005 seems great as a default
    k_regularizer_sz = 0.0005

    if fresh:
        print('  Fresh tissue: using small control point spacing and deformation penalties for photos')
        cp_spacing_2d = 10
        k_regularizer_nonlin = 1.0
    else:
        print('  Fixed tissue: using large control point spacing and deformation penalties for photos')
        cp_spacing_2d = 25.0
        k_regularizer_nonlin = 3.0

    if (not mri_present):
        if (not mesh_present):  # MNI mode
            print('  Using MNI as reference (affine + small control point spacing / low penalty)')
            allow_z_stretch = False
            allow_affine_mri = True
            cp_spacing_3d = 10.0
            k_regularizer_nonlin3d = 1.0

        else:  # mesh mode
            print('  Using mesh as reference')
            allow_z_stretch = True
            allow_affine_mri = None
            cp_spacing_3d = None
            k_regularizer_nonlin3d = 0.0 # irrelevant
    else:
        if (not mesh_present):  # MRI mode
            print('  Using MRI scan as reference (larger control point spacing and penalty)')
            allow_z_stretch = True
            allow_affine_mri = False
            cp_spacing_3d = 20.0
            k_regularizer_nonlin3d = 2.0
        else:  # MRI + mesh mode
            print('  Using mesh as reference and MRI as secondary reference (moderate control point spacing and penalty)')
            allow_z_stretch = True
            allow_affine_mri = True
            cp_spacing_3d = 15.0
            k_regularizer_nonlin3d = 1.5

    if arguments.no_z_stretch:
        allow_z_stretch = False
    
    # override if arugments provided in command line if needed
    arguments.cp_spacing_2d = arguments.cp_spacing_2d if arguments.cp_spacing_2d is not None else cp_spacing_2d
    arguments.k_dif_slice_loss = arguments.k_dif_slice_loss if arguments.k_dif_slice_loss is not None else k_dif_slice_loss
    arguments.k_regularizer = arguments.k_regularizer if arguments.k_regularizer is not None else k_regularizer
    arguments.k_regularizer_nonlin = arguments.k_regularizer_nonlin if arguments.k_regularizer_nonlin is not None else k_regularizer_nonlin
    arguments.k_mesh_loss = arguments.k_mesh_loss if arguments.k_mesh_loss is not None else k_mesh_loss
    arguments.cp_spacing_3d = arguments.cp_spacing_3d if arguments.cp_spacing_3d is not None else cp_spacing_3d
    arguments.k_lncc_mri = arguments.k_lncc_mri if arguments.k_lncc_mri is not None else k_lncc_mri
    arguments.k_dice_mri = arguments.k_dice_mri if arguments.k_dice_mri is not None else k_dice_mri
    arguments.k_regularizer_nonlin3d = arguments.k_regularizer_nonlin3d if arguments.k_regularizer_nonlin3d is not None else k_regularizer_nonlin3d
    arguments.k_regularizer_sz = arguments.k_regularizer_sz if arguments.k_regularizer_sz is not None else k_regularizer_sz


    print('Summary of settings:')
    print('  Control point spacing 2D (in mm):     ' + str(arguments.cp_spacing_2d))
    print('  Weight of slice-to-slice difference:  ' + str(arguments.k_dif_slice_loss))
    print('  Weight of affine regularized (2D):    ' + str(arguments.k_regularizer))
    print('  Weight of nonlinear regularizer (2D): ' + str(arguments.k_regularizer_nonlin))
    if mesh_present:
        print('  Weight of mesh loss:                  ' + str(arguments.k_mesh_loss))
    if mri_present or ((not mri_present) and (not mesh_present)):
        print('  Control point spacing 3D (in mm):     ' + str(arguments.cp_spacing_3d))
        print('  Weight of LNCC with MRI:              ' + str(arguments.k_lncc_mri))
        print('  Weight of Dice with MRI:              ' + str(arguments.k_dice_mri))
        print('  Weight of nonlinear regularizer (3D): ' + str(arguments.k_regularizer_nonlin3d))


    if allow_z_stretch:
        print('  I am allowing adjustment of the slice thickness')
    else:
        print('  I am *not* allowing adjustment of the slice thickness')
    if mri_present:
        if allow_affine_mri:
            print('  I am allowing affine deformation of the MRI')
        else:
            print('  I am *not* allowing affine deformation of the MRI')
    if arguments.weights is None:
        print('  A file with the weights of the slabs was not provided')
    else:
        print('  A file with the weights of the slabs was provided: ' + arguments.weights)
        if (arguments.thickness_cap > 0):
            print('  Maximum allowed thickness (when estimating from weights): ' + str(arguments.thickness_cap))

    return arguments, allow_z_stretch, allow_affine_mri

