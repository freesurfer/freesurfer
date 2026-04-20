from copy import deepcopy

import torch

from batchgeneratorsv2.transforms.base.basic_transform import BasicTransform


class Convert3DTo2DTransform(BasicTransform):
    def apply(self, data_dict, **params):
        if 'image' in data_dict.keys():
            data_dict['nchannels_img'] = deepcopy(data_dict['image']).shape[0]
        if 'segmentation' in data_dict.keys():
            data_dict['nchannels_seg'] = deepcopy(data_dict['segmentation']).shape[0]
        if 'regression_target' in data_dict.keys():
            data_dict['nchannels_regr_trg'] = deepcopy(data_dict['regression_target']).shape[0]
        return super().apply(data_dict, **params)

    def _apply_to_image(self, img: torch.Tensor, **params) -> torch.Tensor:
        shp = img.shape
        return img.reshape((shp[0] * shp[1], *shp[2:]))

    def _apply_to_regr_target(self, regression_target, **params) -> torch.Tensor:
        return self._apply_to_image(regression_target, **params)

    def _apply_to_segmentation(self, segmentation: torch.Tensor, **params) -> torch.Tensor:
        return self._apply_to_image(segmentation, **params)

    def _apply_to_bbox(self, bbox, **params):
        raise NotImplementedError

    def _apply_to_keypoints(self, keypoints, **params):
        raise NotImplementedError


class Convert2DTo3DTransform(BasicTransform):
    def get_parameters(self, **data_dict) -> dict:
        return {i: data_dict[i] for i in
                ['nchannels_img', 'nchannels_seg', 'nchannels_regr_trg']
                if i in data_dict.keys()}

    def apply(self, data_dict, **params):
        data_dict = super().apply(data_dict, **params)
        if 'nchannels_img' in data_dict.keys():
            del data_dict['nchannels_img']
        if 'nchannels_seg' in data_dict.keys():
            del data_dict['nchannels_seg']
        if 'nchannels_regr_trg' in data_dict.keys():
            del data_dict['nchannels_regr_trg']
        return data_dict

    def _apply_to_image(self, img: torch.Tensor, **params) -> torch.Tensor:
        return img.reshape((params['nchannels_img'], img.shape[0] // params['nchannels_img'], *img.shape[1:]))

    def _apply_to_segmentation(self, segmentation: torch.Tensor, **params) -> torch.Tensor:
        return segmentation.reshape(
            (params['nchannels_seg'], segmentation.shape[0] // params['nchannels_seg'], *segmentation.shape[1:]))

    def _apply_to_regr_target(self, regression_target, **params) -> torch.Tensor:
        return regression_target.reshape(
            (params['nchannels_regr_trg'], regression_target.shape[0] // params['nchannels_regr_trg'], *regression_target.shape[1:]))

    def _apply_to_bbox(self, bbox, **params):
        raise NotImplementedError

    def _apply_to_keypoints(self, keypoints, **params):
        raise NotImplementedError


if __name__ == '__main__':
    d = torch.rand((2, 32, 64, 128))
    s = torch.ones((1, 32, 64, 128))

    fwd = Convert3DTo2DTransform()
    bwd = Convert2DTo3DTransform()

    inp = {'image': d, 'segmentation': s}

    tmp = fwd(**inp)
    print(tmp['image'].shape, tmp['segmentation'].shape)
    out = bwd(**tmp)
    print(out['image'].shape, out['segmentation'].shape)
