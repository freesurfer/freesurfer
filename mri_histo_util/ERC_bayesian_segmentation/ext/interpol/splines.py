"""Weights and derivatives of spline orders 0 to 7."""
import torch
from enum import Enum
from .jit_utils import square, cube, pow4, pow5, pow6, pow7


class InterpolationType(Enum):
    nearest = zeroth = 0
    linear = first = 1
    quadratic = second = 2
    cubic = third = 3
    fourth = 4
    fifth = 5
    sixth = 6
    seventh = 7


@torch.jit.script
class Spline:

    def __init__(self, order: int = 1):
        self.order = order

    def weight(self, x):
        w = self.fastweight(x)
        zero = torch.zeros([1], dtype=x.dtype, device=x.device)
        w = torch.where(x.abs() >= (self.order + 1)/2, zero, w)
        return w

    def fastweight(self, x):
        if self.order == 0:
            return torch.ones(x.shape, dtype=x.dtype, device=x.device)
        x = x.abs()
        if self.order == 1:
            return 1 - x
        if self.order == 2:
            x_low = 0.75 - square(x)
            x_up = 0.5 * square(1.5 - x)
            return torch.where(x < 0.5, x_low, x_up)
        if self.order == 3:
            x_low = (x * x * (x - 2.) * 3. + 4.) / 6.
            x_up = cube(2. - x) / 6.
            return torch.where(x < 1., x_low, x_up)
        if self.order == 4:
            x_low = square(x)
            x_low = x_low * (x_low * 0.25 - 0.625) + 115. / 192.
            x_mid = x * (x * (x * (5. - x) / 6. - 1.25) + 5./24.) + 55./96.
            x_up = pow4(x - 2.5) / 24.
            return torch.where(x < 0.5, x_low, torch.where(x < 1.5, x_mid, x_up))
        if self.order == 5:
            x_low = square(x)
            x_low = x_low * (x_low * (0.25 - x / 12.) - 0.5) + 0.55
            x_mid = x * (x * (x * (x * (x / 24. - 0.375) + 1.25) - 1.75) + 0.625) + 0.425
            x_up = pow5(3 - x) / 120.
            return torch.where(x < 1., x_low, torch.where(x < 2., x_mid, x_up))
        if self.order == 6:
            x_low = square(x)
            x_low = x_low * (x_low * (7./48. - x_low/36.) - 77./192.) + 5887./11520.
            x_mid_low = (x * (x * (x * (x * (x * (x / 48. - 7./48.) + 0.328125)
                         - 35./288.) - 91./256.) - 7./768.) + 7861./15360.)
            x_mid_up = (x * (x * (x * (x * (x * (7./60. - x / 120.) - 0.65625)
                        + 133./72.) - 2.5703125) + 1267./960.) + 1379./7680.)
            x_up = pow6(x - 3.5) / 720.
            return torch.where(x < .5, x_low,
                               torch.where(x < 1.5, x_mid_low,
                                           torch.where(x < 2.5, x_mid_up, x_up)))
        if self.order == 7:
            x_low = square(x)
            x_low = (x_low * (x_low * (x_low * (x / 144. - 1./36.)
                     + 1./9.) - 1./3.) + 151./315.)
            x_mid_low = (x * (x * (x * (x * (x * (x * (0.05 - x/240.) - 7./30.)
                         + 0.5) - 7./18.) - 0.1) - 7./90.) + 103./210.)
            x_mid_up = (x * (x * (x * (x * (x * (x * (x / 720. - 1./36.)
                        + 7./30.) - 19./18.) + 49./18.) - 23./6.) + 217./90.)
                        - 139./630.)
            x_up = pow7(4 - x) / 5040.
            return torch.where(x < 1., x_low,
                               torch.where(x < 2., x_mid_low,
                                           torch.where(x < 3., x_mid_up, x_up)))
        raise NotImplementedError

    def grad(self, x):
        if self.order == 0:
            return torch.zeros(x.shape, dtype=x.dtype, device=x.device)
        g = self.fastgrad(x)
        zero = torch.zeros([1], dtype=x.dtype, device=x.device)
        g = torch.where(x.abs() >= (self.order + 1)/2, zero, g)
        return g

    def fastgrad(self, x):
        if self.order == 0:
            return torch.zeros(x.shape, dtype=x.dtype, device=x.device)
        return self._fastgrad(x.abs()).mul(x.sign())

    def _fastgrad(self, x):
        if self.order == 1:
            return torch.ones(x.shape, dtype=x.dtype, device=x.device)
        if self.order == 2:
            return torch.where(x < 0.5, -2*x, x - 1.5)
        if self.order == 3:
            g_low = x * (x * 1.5 - 2)
            g_up = -0.5 * square(2 - x)
            return torch.where(x < 1, g_low, g_up)
        if self.order == 4:
            g_low = x * (square(x) - 1.25)
            g_mid = x * (x * (x * (-2./3.) + 2.5) - 2.5) + 5./24.
            g_up = cube(2. * x - 5.) / 48.
            return torch.where(x < 0.5, g_low,
                               torch.where(x < 1.5, g_mid, g_up))
        if self.order == 5:
            g_low = x * (x * (x * (x * (-5./12.) + 1.)) - 1.)
            g_mid = x * (x * (x * (x * (5./24.) - 1.5) + 3.75) - 3.5) + 0.625
            g_up = pow4(x - 3.) / (-24.)
            return torch.where(x < 1, g_low,
                               torch.where(x < 2, g_mid, g_up))
        if self.order == 6:
            g_low = square(x)
            g_low = x * (g_low * (7./12.) - square(g_low) / 6. - 77./96.)
            g_mid_low = (x * (x * (x * (x * (x * 0.125 - 35./48.) + 1.3125)
                         - 35./96.) - 0.7109375) - 7./768.)
            g_mid_up = (x * (x * (x * (x * (x / (-20.) + 7./12.) - 2.625)
                        + 133./24.) - 5.140625) + 1267./960.)
            g_up = pow5(2*x - 7) / 3840.
            return torch.where(x < 0.5, g_low,
                               torch.where(x < 1.5, g_mid_low,
                                           torch.where(x < 2.5, g_mid_up,
                                                       g_up)))
        if self.order == 7:
            g_low = square(x)
            g_low = x * (g_low * (g_low * (x * (7./144.) - 1./6.) + 4./9.) - 2./3.)
            g_mid_low = (x * (x * (x * (x * (x * (x * (-7./240.) + 3./10.)
                         - 7./6.) + 2.) - 7./6.) - 1./5.) - 7./90.)
            g_mid_up = (x * (x * (x * (x * (x * (x * (7./720.) - 1./6.)
                        + 7./6.) - 38./9.) + 49./6.) - 23./3.) + 217./90.)
            g_up = pow6(x - 4) / (-720.)
            return torch.where(x < 1, g_low,
                               torch.where(x < 2, g_mid_low,
                                           torch.where(x < 3, g_mid_up, g_up)))
        raise NotImplementedError

    def hess(self, x):
        if self.order == 0:
            return torch.zeros(x.shape, dtype=x.dtype, device=x.device)
        h = self.fasthess(x)
        zero = torch.zeros([1], dtype=x.dtype, device=x.device)
        h = torch.where(x.abs() >= (self.order + 1)/2, zero, h)
        return h

    def fasthess(self, x):
        if self.order in (0, 1):
            return torch.zeros(x.shape, dtype=x.dtype, device=x.device)
        x = x.abs()
        if self.order == 2:
            one = torch.ones([1], dtype=x.dtype, device=x.device)
            return torch.where(x < 0.5, -2 * one, one)
        if self.order == 3:
            return torch.where(x < 1, 3. * x - 2., 2. - x)
        if self.order == 4:
            return torch.where(x < 0.5, 3. * square(x) - 1.25,
                               torch.where(x < 1.5, x * (-2. * x + 5.) - 2.5,
                                           square(2. * x - 5.) / 8.))
        if self.order == 5:
            h_low = square(x)
            h_low = - h_low * (x * (5./3.) - 3.) - 1.
            h_mid = x * (x * (x * (5./6.) - 9./2.) + 15./2.) - 7./2.
            h_up = 9./2. - x * (x * (x/6. - 3./2.) + 9./2.)
            return torch.where(x < 1, h_low,
                               torch.where(x < 2, h_mid, h_up))
        if self.order == 6:
            h_low = square(x)
            h_low = - h_low * (h_low * (5./6) - 7./4.) - 77./96.
            h_mid_low = (x * (x * (x * (x * (5./8.) - 35./12.) + 63./16.)
                         - 35./48.) - 91./128.)
            h_mid_up = -(x * (x * (x * (x/4. - 7./3.) + 63./8.) - 133./12.)
                         + 329./64.)
            h_up = (x * (x * (x * (x/24. - 7./12.) + 49./16.) - 343./48.)
                    + 2401./384.)
            return torch.where(x < 0.5, h_low,
                               torch.where(x < 1.5, h_mid_low,
                                           torch.where(x < 2.5, h_mid_up,
                                                       h_up)))
        if self.order == 7:
            h_low = square(x)
            h_low = h_low * (h_low*(x * (7./24.) - 5./6.) + 4./3.) - 2./3.
            h_mid_low = - (x * (x * (x * (x * (x * (7./40.) - 3./2.) + 14./3.)
                           - 6.) + 7./3.) + 1./5.)
            h_mid_up = (x * (x * (x * (x * (x * (7./120.) - 5./6.) + 14./3.)
                        - 38./3.) + 49./3.) - 23./3.)
            h_up = - (x * (x * (x * (x * (x/120. - 1./6.) + 4./3.) - 16./3.)
                      + 32./3.) - 128./15.)
            return torch.where(x < 1, h_low,
                               torch.where(x < 2, h_mid_low,
                                           torch.where(x < 3, h_mid_up,
                                                       h_up)))
        raise NotImplementedError

