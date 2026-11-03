"""Independent, unilateral distributed contact model for the contact-only pilot.

This module deliberately imports no benchmark, footprint_kernels or scene generator.
Pressure is proportional to A(x+xi) * [d-xi**2/(2R)]_+**n. Its positive elastic
potential gives force and pressure moment by integration. p=n+1/2 is the uniform
force exponent; p=1 is a semielliptic cylindrical approximation and p=1.5 a
parabolic elastic-foundation approximation. Neither is an exact heterogeneous
half-space solution. The fitted scale has units N/m**p, NOT Young's modulus.

The probe radius is a declared geometric assumption (28 mm), not a plant read.
Two smooth edges describe one local patch; they do not assume scene ridge locations.
"""
from dataclasses import dataclass
import numpy as np
from scipy.special import roots_jacobi


@dataclass(frozen=True)
class Geometry:
    radius: float = .028
    exponent: float = 1.0
    edge: float = .0013
    quadrature: int = 32


class FinitePatch:
    def __init__(self, geometry=Geometry(), spatial=False):
        if not .5 < geometry.exponent <= 2 or geometry.radius <= 0:
            raise ValueError('positive radius and integrable elastic pressure required')
        self.g, self.spatial = geometry, bool(spatial)
        n = geometry.exponent-.5
        self.u, self.w = roots_jacobi(geometry.quadrature, n, n)
        self.w /= self.w.sum()
        self.uk, self.wk = roots_jacobi(geometry.quadrature, n-1, n-1)
        # Integral n*(1-u^2)^(n-1) / integral (1-u^2)^n equals p.
        self.wk *= n / roots_jacobi(geometry.quadrature, n, n)[1].sum()

    def material(self, xx, theta):
        if not self.spatial:
            return np.ones_like(xx)
        _, lr, b, width = theta
        left = 1/(1+np.exp(np.clip(-(xx-b)/self.g.edge, -50, 50)))
        right = 1/(1+np.exp(np.clip((xx-b-width)/self.g.edge, -50, 50)))
        # Positive interpolation between two positive foundation coefficients.
        return 1+np.expm1(lr)*left*right

    def response(self, x, d, theta, stiffness=True):
        x, d = np.broadcast_arrays(np.asarray(x, float), np.asarray(d, float))
        d = np.maximum(d, 0.)
        a = np.sqrt(2*self.g.radius*d)
        A = self.material(x[..., None]+a[..., None]*self.u, theta)
        c = np.exp(theta[0]); scale = c*d**self.g.exponent
        f = scale*(A@self.w)
        m = scale*a*(A@(self.w*self.u))
        k = None
        if stiffness:
            Ak = self.material(x[..., None]+a[..., None]*self.uk, theta)
            k = c*np.maximum(d, 1e-12)**(self.g.exponent-1)*(Ak@self.wk)
            k = np.where(d > 0, k, 0.)
        return f, m, k

    def inverse(self, x, f, theta):
        """Solve F(d)=observed load using positive derivative; no measured depth."""
        f = np.maximum(np.asarray(f, float), .005)
        d = np.clip((f/np.exp(theta[0]))**(1/self.g.exponent), 1e-7, .06)
        for _ in range(16):
            fe, _, k = self.response(x, d, theta)
            if np.all(np.abs(fe-f)<=1e-13*np.maximum(f,1.)):
                break
            step = (fe-f)/np.maximum(k, 1.)
            d = np.clip(d-np.clip(step, -.6*d, .6*d), 1e-8, .08)
        return d


def orthogonal_residual(B, y):
    """Rank-aware elimination, including rank-deficient wrench nuisance columns."""
    U, s, _ = np.linalg.svd(B, full_matrices=False)
    keep = s > max(B.shape)*np.finfo(float).eps*max(s[0], 1.)
    Q = U[:, keep]
    return y-Q@(Q.T@y)
