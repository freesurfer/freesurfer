import torch
from torch import nn

from dynamic_network_architectures.building_blocks.plain_conv_encoder import PlainConvEncoder
from dynamic_network_architectures.building_blocks.helper import get_matching_pool_op, get_default_network_config

_VGG_CONFIGS = {
    '16': {'features_per_stage': (64, 128, 256, 512, 512, 512), 'n_conv_per_stage': (2, 2, 2, 3, 3, 3),
           'strides': (1, 2, 2, 2, 2, 2)},
    '19': {'features_per_stage': (64, 128, 256, 512, 512, 512), 'n_conv_per_stage': (2, 2, 3, 3, 4, 4),
           'strides': (1, 2, 2, 2, 2, 2)},
    '16_cifar': {'features_per_stage': (64, 128, 256, 512), 'n_conv_per_stage': (2, 3, 5, 5), 'strides': (1, 2, 2, 2)},
    '19_cifar': {'features_per_stage': (64, 128, 256, 512), 'n_conv_per_stage': (3, 4, 5, 6), 'strides': (1, 2, 2, 2)},
}

_VGG_OPS = {
    1: {'conv_op': nn.Conv1d, 'norm_op': nn.BatchNorm1d},
    2: {'conv_op': nn.Conv2d, 'norm_op': nn.BatchNorm2d},
    3: {'conv_op': nn.Conv3d, 'norm_op': nn.BatchNorm3d},
}


class VGG(nn.Module):
    def __init__(self, n_classes: int, n_input_channel: int = 3, config='16', input_dimension=2):
        """
        This is not 1:1 VGG because it does not have the bloated fully connected layers at the end. Since these were
        counted towards the XX layers as well, we increase the number of convolutional layers so that we have the
        desired number of conv layers in total

        We also use batchnorm
        """
        super().__init__()
        cfg = _VGG_CONFIGS[config]
        ops = get_default_network_config(dimension=input_dimension)
        self.encoder = PlainConvEncoder(
            n_input_channel, n_stages=len(cfg['features_per_stage']), features_per_stage=cfg['features_per_stage'],
            conv_op=ops['conv_op'],
            kernel_sizes=3, strides=cfg['strides'], n_conv_per_stage=cfg['n_conv_per_stage'], conv_bias=False,
            norm_op=ops['norm_op'], norm_op_kwargs=None, dropout_op=None, dropout_op_kwargs=None, nonlin=nn.ReLU,
            nonlin_kwargs={'inplace': True}, return_skips=False
        )
        self.gap = get_matching_pool_op(conv_op=ops['conv_op'], adaptive=True, pool_type='avg')(1)
        self.classifier = nn.Linear(cfg['features_per_stage'][-1], n_classes, True)

    def forward(self, x):
        x = self.encoder(x)
        x = self.gap(x).squeeze()
        return self.classifier(x)

    def compute_conv_feature_map_size(self, input_size):
        return self.encoder.compute_conv_feature_map_size(input_size)


class VGG16(VGG):
    def __init__(self, n_classes: int, n_input_channel: int = 3, input_dimension: int = 2):
        super().__init__(n_classes, n_input_channel, config='16', input_dimension=input_dimension)


class VGG19(VGG):
    def __init__(self, n_classes: int, n_input_channel: int = 3, input_dimension: int = 2):
        super().__init__(n_classes, n_input_channel, config='19', input_dimension=input_dimension)


class VGG16_cifar(VGG):
    def __init__(self, n_classes: int, n_input_channel: int = 3, input_dimension: int = 2):
        super().__init__(n_classes, n_input_channel, config='16_cifar', input_dimension=input_dimension)


class VGG19_cifar(VGG):
    def __init__(self, n_classes: int, n_input_channel: int = 3, input_dimension: int = 2):
        super().__init__(n_classes, n_input_channel, config='19_cifar', input_dimension=input_dimension)


if __name__ == '__main__':
    data = torch.rand((1, 3, 32, 32))

    model = VGG19_cifar(10, 3)
    import hiddenlayer as hl

    g = hl.build_graph(model, data,
                       transforms=None)
    g.save("network_architecture.pdf")
    del g

    print(model.compute_conv_feature_map_size((32, 32)))