import torch
from dynamic_network_architectures.building_blocks.residual_encoders import ResidualEncoder, BottleneckD, BasicBlockD
from dynamic_network_architectures.building_blocks.helper import get_matching_pool_op, get_default_network_config
from dynamic_network_architectures.building_blocks.simple_conv_blocks import ConvDropoutNormReLU
from torch import nn

_ResNet_CONFIGS = {
    '18': {'features_per_stage': (64, 128, 256, 512), 'n_blocks_per_stage': (2, 2, 2, 2), 'strides': (1, 2, 2, 2),
           'block': BasicBlockD, 'bottleneck_channels': None, 'disable_default_stem': True, 'stem_channels': None},
    '34': {'features_per_stage': (64, 128, 256, 512), 'n_blocks_per_stage': (3, 4, 6, 3), 'strides': (1, 2, 2, 2),
           'block': BasicBlockD, 'bottleneck_channels': None, 'disable_default_stem': True, 'stem_channels': None},
    '50': {'features_per_stage': (64, 128, 256, 512), 'n_blocks_per_stage': (4, 6, 10, 5), 'strides': (1, 2, 2, 2),
           'block': BasicBlockD, 'bottleneck_channels': None, 'disable_default_stem': True, 'stem_channels': None},
    '152': {'features_per_stage': (64, 128, 256, 512), 'n_blocks_per_stage': (4, 13, 55, 4), 'strides': (1, 2, 2, 2),
            'block': BasicBlockD, 'bottleneck_channels': None, 'disable_default_stem': True, 'stem_channels': None},
    '50_bn': {'features_per_stage': (256, 512, 1024, 2048), 'n_blocks_per_stage': (3, 4, 6, 3), 'strides': (1, 2, 2, 2),
              'block': BottleneckD, 'bottleneck_channels': (64, 128, 256, 512), 'disable_default_stem': True,
              'stem_channels': 64},
    '152_bn': {'features_per_stage': (256, 512, 1024, 2048), 'n_blocks_per_stage': (3, 8, 36, 3),
               'strides': (1, 2, 2, 2),
               'block': BottleneckD, 'bottleneck_channels': (64, 128, 256, 512), 'disable_default_stem': True,
               'stem_channels': 64},
    '18_cifar': {'features_per_stage': (64, 128, 256, 512), 'n_blocks_per_stage': (2, 2, 2, 2), 'strides': (1, 2, 2, 2),
                 'block': BasicBlockD, 'bottleneck_channels': None, 'disable_default_stem': False,
                 'stem_channels': None},
    '34_cifar': {'features_per_stage': (64, 128, 256, 512), 'n_blocks_per_stage': (3, 4, 6, 3), 'strides': (1, 2, 2, 2),
                 'block': BasicBlockD, 'bottleneck_channels': None, 'disable_default_stem': False,
                 'stem_channels': None},
    '50_cifar': {'features_per_stage': (64, 128, 256, 512), 'n_blocks_per_stage': (4, 6, 10, 5),
                 'strides': (1, 2, 2, 2),
                 'block': BasicBlockD, 'bottleneck_channels': None, 'disable_default_stem': False,
                 'stem_channels': None},
    '152_cifar': {'features_per_stage': (64, 128, 256, 512), 'n_blocks_per_stage': (4, 13, 55, 4),
                  'strides': (1, 2, 2, 2),
                  'block': BasicBlockD, 'bottleneck_channels': None, 'disable_default_stem': False,
                  'stem_channels': None},
    '50_cifar_bn': {'features_per_stage': (256, 512, 1024, 2048), 'n_blocks_per_stage': (3, 4, 6, 3),
                    'strides': (1, 2, 2, 2),
                    'block': BottleneckD, 'bottleneck_channels': (64, 128, 256, 512), 'disable_default_stem': False,
                    'stem_channels': 64},
    '152_cifar_bn': {'features_per_stage': (256, 512, 1024, 2048), 'n_blocks_per_stage': (3, 8, 36, 3),
                     'strides': (1, 2, 2, 2),
                     'block': BottleneckD, 'bottleneck_channels': (64, 128, 256, 512), 'disable_default_stem': False,
                     'stem_channels': 64},
}


class ResNetD(nn.Module):
    def __init__(self, n_classes: int, n_input_channel: int = 3, config='18', input_dimension=2,
                 final_layer_dropout=0.0, stochastic_depth_p=0.0, squeeze_excitation=False,
                 squeeze_excitation_rd_ratio=1./16):
        """
        Implements ResNetD (https://arxiv.org/pdf/1812.01187.pdf).
        Args:
            n_classes: Number of classes
            n_input_channel: Number of input channels (e.g. 3 for RGB)
            config: Configuration of the ResNet
            input_dimension: Number of dimensions of the data (1, 2 or 3)
            final_layer_dropout: Probability of dropout before the final classifier
            stochastic_depth_p: Stochastic Depth probability
            squeeze_excitation: Whether Squeeze and Excitation should be applied
            squeeze_excitation_rd_ratio: Squeeze and Excitation Reduction Ratio
        Returns:
            ResNet Model
        """
        super().__init__()
        self.input_channels = n_input_channel
        self.cfg = _ResNet_CONFIGS[config]
        self.ops = get_default_network_config(dimension=input_dimension)
        self.final_layer_dropout_p = final_layer_dropout

        if self.cfg['disable_default_stem']:
            stem_features = self.cfg['stem_channels'] if self.cfg['stem_channels'] is not None else \
            self.cfg['features_per_stage'][0]
            self.stem = self._build_imagenet_stem_D(stem_features)
            encoder_input_features = stem_features
        else:
            encoder_input_features = n_input_channel
            self.stem = None

        self.encoder = ResidualEncoder(encoder_input_features, n_stages=len(self.cfg['features_per_stage']),
                                       features_per_stage=self.cfg['features_per_stage'], conv_op=self.ops['conv_op'],
                                       kernel_sizes=3, strides=self.cfg['strides'],
                                       n_blocks_per_stage=self.cfg['n_blocks_per_stage'], conv_bias=False,
                                       norm_op=self.ops['norm_op'], norm_op_kwargs=None, dropout_op=None,
                                       dropout_op_kwargs=None, nonlin=nn.ReLU,
                                       nonlin_kwargs={'inplace': True}, block=self.cfg['block'],
                                       bottleneck_channels=self.cfg['bottleneck_channels'], return_skips=False,
                                       disable_default_stem=self.cfg['disable_default_stem'],
                                       stem_channels=self.cfg['stem_channels'],
                                       stochastic_depth_p=stochastic_depth_p,
                                       squeeze_excitation=squeeze_excitation,
                                       squeeze_excitation_reduction_ratio=squeeze_excitation_rd_ratio)

        self.gap = get_matching_pool_op(conv_op=self.ops['conv_op'], adaptive=True, pool_type='avg')(1)
        self.classifier = nn.Linear(self.cfg['features_per_stage'][-1], n_classes, True)
        self.final_layer_dropout = self.ops['dropout_op'](p=self.final_layer_dropout_p)

    def forward(self, x):
        if self.stem is not None:
            x = self.stem(x)
        x = self.encoder(x)
        x = self.gap(x)
        x = self.final_layer_dropout(x).squeeze()

        return self.classifier(x)

    def _build_imagenet_stem_D(self, stem_features):
        """
        https://arxiv.org/pdf/1812.01187.pdf

        use 3 3x3(x3) convs instead of one 7x7. Stride is located in first conv.

        Fig2 b) describes this
        :return:
        """
        c1 = ConvDropoutNormReLU(self.ops['conv_op'], self.input_channels, stem_features, 3, 2, False,
                                 self.ops['norm_op'], None, None, None, nn.ReLU, {'inplace': True})
        c2 = ConvDropoutNormReLU(self.ops['conv_op'], stem_features, stem_features, 3, 1, False,
                                 self.ops['norm_op'], None, None, None, nn.ReLU, {'inplace': True})
        c3 = ConvDropoutNormReLU(self.ops['conv_op'], stem_features, stem_features, 3, 1, False,
                                 self.ops['norm_op'], None, None, None, nn.ReLU, {'inplace': True})
        pl = get_matching_pool_op(conv_op=self.ops['conv_op'], adaptive=False, pool_type='max')(2)
        stem = nn.Sequential(c1, c2, c3, pl)
        return stem


class ResNet18_CIFAR(ResNetD):
    def __init__(self, n_classes: int, n_input_channels: int = 3, input_dimension: int = 2,
                 final_layer_dropout: float = 0.0, stochastic_depth_p: float = 0.0, squeeze_excitation: bool = False,
                 squeeze_excitation_rd_ratio: float = 1./16):
        super().__init__(n_classes, n_input_channels, config='18_cifar', input_dimension=input_dimension,
                         final_layer_dropout=final_layer_dropout, stochastic_depth_p=stochastic_depth_p,
                         squeeze_excitation=squeeze_excitation, squeeze_excitation_rd_ratio=squeeze_excitation_rd_ratio)

class ResNet34_CIFAR(ResNetD):
    def __init__(self, n_classes: int, n_input_channels: int = 3, input_dimension: int = 2,
                 final_layer_dropout: float = 0.0, stochastic_depth_p: float = 0.0, squeeze_excitation: bool = False,
                 squeeze_excitation_rd_ratio: float = 1./16):
        super().__init__(n_classes, n_input_channels, config='34_cifar', input_dimension=input_dimension,
                         final_layer_dropout=final_layer_dropout, stochastic_depth_p=stochastic_depth_p,
                         squeeze_excitation=squeeze_excitation, squeeze_excitation_rd_ratio=squeeze_excitation_rd_ratio)

class ResNet50_CIFAR(ResNetD):
    def __init__(self, n_classes: int, n_input_channels: int = 3, input_dimension: int = 2,
                 final_layer_dropout: float = 0.0, stochastic_depth_p: float = 0.0, squeeze_excitation: bool = False,
                 squeeze_excitation_rd_ratio: float = 1./16):
        super().__init__(n_classes, n_input_channels, config='50_cifar', input_dimension=input_dimension,
                         final_layer_dropout=final_layer_dropout, stochastic_depth_p=stochastic_depth_p,
                         squeeze_excitation=squeeze_excitation, squeeze_excitation_rd_ratio=squeeze_excitation_rd_ratio)

class ResNet152_CIFAR(ResNetD):
    def __init__(self, n_classes: int, n_input_channels: int = 3, input_dimension: int = 2,
                 final_layer_dropout: float = 0.0, stochastic_depth_p: float = 0.0, squeeze_excitation: bool = False,
                 squeeze_excitation_rd_ratio: float = 1./16):
        super().__init__(n_classes, n_input_channels, config='152_cifar', input_dimension=input_dimension,
                         final_layer_dropout=final_layer_dropout, stochastic_depth_p=stochastic_depth_p,
                         squeeze_excitation=squeeze_excitation, squeeze_excitation_rd_ratio=squeeze_excitation_rd_ratio)

class ResNet50bn_CIFAR(ResNetD):
    def __init__(self, n_classes: int, n_input_channels: int = 3, input_dimension: int = 2,
                 final_layer_dropout: float = 0.0, stochastic_depth_p: float = 0.0, squeeze_excitation: bool = False,
                 squeeze_excitation_rd_ratio: float = 1./16):
        super().__init__(n_classes, n_input_channels, config='50_cifar_bn', input_dimension=input_dimension,
                         final_layer_dropout=final_layer_dropout, stochastic_depth_p=stochastic_depth_p,
                         squeeze_excitation=squeeze_excitation, squeeze_excitation_rd_ratio=squeeze_excitation_rd_ratio)

class ResNet152bn_CIFAR(ResNetD):
    def __init__(self, n_classes: int, n_input_channels: int = 3, input_dimension: int = 2,
                 final_layer_dropout: float = 0.0, stochastic_depth_p: float = 0.0, squeeze_excitation: bool = False,
                 squeeze_excitation_rd_ratio: float = 1./16):
        super().__init__(n_classes, n_input_channels, config='152_cifar_bn', input_dimension=input_dimension,
                         final_layer_dropout=final_layer_dropout, stochastic_depth_p=stochastic_depth_p,
                         squeeze_excitation=squeeze_excitation, squeeze_excitation_rd_ratio=squeeze_excitation_rd_ratio)

class ResNet18(ResNetD):
    def __init__(self, n_classes: int, n_input_channels: int = 3, input_dimension: int = 2,
                 final_layer_dropout: float = 0.0, stochastic_depth_p: float = 0.0, squeeze_excitation: bool = False,
                 squeeze_excitation_rd_ratio: float = 1./16):
        super().__init__(n_classes, n_input_channels, config='18', input_dimension=input_dimension,
                         final_layer_dropout=final_layer_dropout, stochastic_depth_p=stochastic_depth_p,
                         squeeze_excitation=squeeze_excitation, squeeze_excitation_rd_ratio=squeeze_excitation_rd_ratio)

class ResNet34(ResNetD):
    def __init__(self, n_classes: int, n_input_channels: int = 3, input_dimension: int = 2,
                 final_layer_dropout: float = 0.0, stochastic_depth_p: float = 0.0, squeeze_excitation: bool = False,
                 squeeze_excitation_rd_ratio: float = 1./16):
        super().__init__(n_classes, n_input_channels, config='34', input_dimension=input_dimension,
                         final_layer_dropout=final_layer_dropout, stochastic_depth_p=stochastic_depth_p,
                         squeeze_excitation=squeeze_excitation, squeeze_excitation_rd_ratio=squeeze_excitation_rd_ratio)

class ResNet50(ResNetD):
    def __init__(self, n_classes: int, n_input_channels: int = 3, input_dimension: int = 2,
                 final_layer_dropout: float = 0.0, stochastic_depth_p: float = 0.0, squeeze_excitation: bool = False,
                 squeeze_excitation_rd_ratio: float = 1./16):
        super().__init__(n_classes, n_input_channels, config='50', input_dimension=input_dimension,
                         final_layer_dropout=final_layer_dropout, stochastic_depth_p=stochastic_depth_p,
                         squeeze_excitation=squeeze_excitation, squeeze_excitation_rd_ratio=squeeze_excitation_rd_ratio)

class ResNet152(ResNetD):
    def __init__(self, n_classes: int, n_input_channels: int = 3, input_dimension: int = 2,
                 final_layer_dropout: float = 0.0, stochastic_depth_p: float = 0.0, squeeze_excitation: bool = False,
                 squeeze_excitation_rd_ratio: float = 1./16):
        super().__init__(n_classes, n_input_channels, config='152', input_dimension=input_dimension,
                         final_layer_dropout=final_layer_dropout, stochastic_depth_p=stochastic_depth_p,
                         squeeze_excitation=squeeze_excitation, squeeze_excitation_rd_ratio=squeeze_excitation_rd_ratio)

class ResNet50bn(ResNetD):
    def __init__(self, n_classes: int, n_input_channels: int = 3, input_dimension: int = 2,
                 final_layer_dropout: float = 0.0, stochastic_depth_p: float = 0.0, squeeze_excitation: bool = False,
                 squeeze_excitation_rd_ratio: float = 1./16):
        super().__init__(n_classes, n_input_channels, config='50_bn', input_dimension=input_dimension,
                         final_layer_dropout=final_layer_dropout, stochastic_depth_p=stochastic_depth_p,
                         squeeze_excitation=squeeze_excitation, squeeze_excitation_rd_ratio=squeeze_excitation_rd_ratio)

class ResNet152bn(ResNetD):
    def __init__(self, n_classes: int, n_input_channels: int = 3, input_dimension: int = 2,
                 final_layer_dropout: float = 0.0, stochastic_depth_p: float = 0.0, squeeze_excitation: bool = False,
                 squeeze_excitation_rd_ratio: float = 1./16):
        super().__init__(n_classes, n_input_channels, config='152_bn', input_dimension=input_dimension,
                         final_layer_dropout=final_layer_dropout, stochastic_depth_p=stochastic_depth_p,
                         squeeze_excitation=squeeze_excitation, squeeze_excitation_rd_ratio=squeeze_excitation_rd_ratio)


if __name__ == '__main__':
    data = torch.rand((1, 3, 224, 224))

    model = ResNet50bn(10, 3)
    import hiddenlayer as hl

    g = hl.build_graph(model, data,
                       transforms=None)
    g.save("network_architecture.pdf")
    del g

    #print(model.compute_conv_feature_map_size((32, 32)))