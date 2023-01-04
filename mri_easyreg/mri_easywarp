import os
import argparse
import numpy as np
import torch
import surfa as sf
import nibabel as nib


def main():

    parser = argparse.ArgumentParser(description="EasyReg (warping code): deep learning registration simple and easy", epilog='\n')

    # input/outputs
    parser.add_argument("--i", help="Input image")
    parser.add_argument("--o", help="Output (deformed) image")
    parser.add_argument("--field", help="Deformation field")
    parser.add_argument("--nearest", action="store_true", help="(optional) Use nearest neighbor (rather than linear) interpolation")
    parser.add_argument("--threads", type=int, default=1, help="(optional) Number of cores to be used. Default is 1. You can use -1 to use all available cores")

    # parse commandline
    args = parser.parse_args()

    #############

    if args.i is None:
        sf.system.fatal('Input image must be provided')
    if args.o is None:
        sf.system.fatal('Output (deformed) image must be provided')
    if args.field is None:
        sf.system.fatal('Deformation field must be provided')

    # limit the number of threads to be used if running on CPU
    os.environ['CUDA_VISIBLE_DEVICES'] = '-1'
    if args.threads == 1:
        print('using 1 thread')
    elif args.threads<0:
        args.threads = os.cpu_count()
        print('using all available threads ( %s )' % args.threads)
    else:
        print('using %s threads' % args.threads)
    torch.set_num_threads(args.threads)

    print('Reading deformation field')
    field_buffer, field_aff, field_h = load_volume(args.field, im_only=False, squeeze=True, dtype=None)
    if len(field_buffer.shape) !=4:
        sf.system.fatal('field must be 4D array')
    if field_buffer.shape[3] != 3:
        sf.system.fatal('field must have 3 frames')

    print('Reading input image')
    input_buffer, input_aff, input_h = load_volume(args.i, im_only=False, squeeze=True, dtype=None)

    print('Deforming (interpolating)')
    affine = torch.tensor(np.linalg.inv(input_aff), device='cpu')
    field_buffer = torch.tensor(field_buffer, device='cpu')
    II = affine[0, 0] * field_buffer[:,:,:,0]  + affine[0, 1] * field_buffer[:,:,:,1]  + affine[0, 2] * field_buffer[:,:,:,2]  + affine[0, 3]
    JJ = affine[1, 0] * field_buffer[:,:,:,0]  + affine[1, 1] * field_buffer[:,:,:,1]  + affine[1, 2] * field_buffer[:,:,:,2]  + affine[1, 3]
    KK = affine[2, 0] * field_buffer[:,:,:,0]  + affine[2, 1] * field_buffer[:,:,:,1]  + affine[2, 2] * field_buffer[:,:,:,2]  + affine[2, 3]

    if args.nearest:
        Y = fast_3D_interp_torch(torch.tensor(input_buffer, device='cpu', requires_grad=False), II, JJ, KK, 'nearest')
    else:
        Y = fast_3D_interp_torch(torch.tensor(input_buffer, device='cpu', requires_grad=False), II, JJ, KK, 'linear')

    print('Saving to disk')
    save_volume(Y.numpy(), field_aff, field_h, args.o)

    print('All done!')



#######################
# Auxiliary functions #
#######################

def load_volume(path_volume, im_only=True, squeeze=True, dtype=None):

    assert path_volume.endswith(('.nii', '.nii.gz', '.mgz', '.npz')), 'Unknown data file: %s' % path_volume

    if path_volume.endswith(('.nii', '.nii.gz', '.mgz')):
        x = nib.load(path_volume)
        if squeeze:
            volume = np.squeeze(x.get_fdata())
        else:
            volume = x.get_fdata()
        aff = x.affine
        header = x.header
    else:  # npz
        volume = np.load(path_volume)['vol_data']
        if squeeze:
            volume = np.squeeze(volume)
        aff = np.eye(4)
        header = nib.Nifti1Header()
    if dtype is not None:
        if 'int' in dtype:
            volume = np.round(volume)
        volume = volume.astype(dtype=dtype)

    if im_only:
        return volume
    else:
        return volume, aff, header


def save_volume(volume, aff, header, path):
    mkdir(os.path.dirname(path))
    if '.npz' in path:
        np.savez_compressed(path, vol_data=volume)
    else:
        if header is None:
            header = nib.Nifti1Header()
        if isinstance(aff, str):
            if aff == 'FS':
                aff = np.array([[-1, 0, 0, 0], [0, 0, 1, 0], [0, -1, 0, 0], [0, 0, 0, 1]])
        elif aff is None:
            aff = np.eye(4)
        nifty = nib.Nifti1Image(volume, aff, header)

        nib.save(nifty, path)


def mkdir(path_dir):

    if path_dir[-1] == '/':
        path_dir = path_dir[:-1]
    if not os.path.isdir(path_dir):
        list_dir_to_create = [path_dir]
        while not os.path.isdir(os.path.dirname(list_dir_to_create[-1])):
            list_dir_to_create.append(os.path.dirname(list_dir_to_create[-1]))
        for dir_to_create in reversed(list_dir_to_create):
            os.mkdir(dir_to_create)



def fast_3D_interp_torch(X, II, JJ, KK, mode):
    if mode=='nearest':
        IIr = torch.round(II).long()
        JJr = torch.round(JJ).long()
        KKr = torch.round(KK).long()
        IIr[IIr < 0] = 0
        JJr[JJr < 0] = 0
        KKr[KKr < 0] = 0
        IIr[IIr > (X.shape[0] - 1)] = (X.shape[0] - 1)
        JJr[JJr > (X.shape[1] - 1)] = (X.shape[1] - 1)
        KKr[KKr > (X.shape[2] - 1)] = (X.shape[2] - 1)
        Y = X[IIr, JJr, KKr]
    elif mode=='linear':
        ok = (II>0) & (JJ>0) & (KK>0) & (II<=X.shape[0]-1) & (JJ<=X.shape[1]-1) & (KK<=X.shape[2]-1)
        IIv = II[ok]
        JJv = JJ[ok]
        KKv = KK[ok]

        fx = torch.floor(IIv).long()
        cx = fx + 1
        cx[cx > (X.shape[0] - 1)] = (X.shape[0] - 1)
        wcx = IIv - fx
        wfx = 1 - wcx

        fy = torch.floor(JJv).long()
        cy = fy + 1
        cy[cy > (X.shape[1] - 1)] = (X.shape[1] - 1)
        wcy = JJv - fy
        wfy = 1 - wcy

        fz = torch.floor(KKv).long()
        cz = fz + 1
        cz[cz > (X.shape[2] - 1)] = (X.shape[2] - 1)
        wcz = KKv - fz
        wfz = 1 - wcz

        c000 = X[fx, fy, fz]
        c100 = X[cx, fy, fz]
        c010 = X[fx, cy, fz]
        c110 = X[cx, cy, fz]
        c001 = X[fx, fy, cz]
        c101 = X[cx, fy, cz]
        c011 = X[fx, cy, cz]
        c111 = X[cx, cy, cz]

        c00 = c000 * wfx + c100 * wcx
        c01 = c001 * wfx + c101 * wcx
        c10 = c010 * wfx + c110 * wcx
        c11 = c011 * wfx + c111 * wcx

        c0 = c00 * wfy + c10 * wcy
        c1 = c01 * wfy + c11 * wcy

        c = c0 * wfz + c1 * wcz

        Y = torch.zeros(II.shape, device='cpu')
        Y[ok] = c.float()

    else:
        sf.system.fatal('mode must be linear or nearest')

    return Y




# execute script
if __name__ == '__main__':
    main()