"""E3/B0-related cylindrical pressure candidate in the SAME inverse-depth interface.

Independent probe geometry R=28 mm, outer displacement reference H=25 mm;
no benchmark imports. Uses the E3 reference-depth approximation, not the plant's
psi function. F = C a^2 <A>, d=a^2/(2R)[ln(4H/a)+1/2]. Positive pressure,
positive stiffness; C absorbs probe length and elastic coefficient. This retains
finite-footprint physics while changing the E3 global single-edge integration.
"""
import math
import numpy as np
from finite_patch_model import FinitePatch, Geometry


class LinePatch(FinitePatch):
    reference = .025

    def __init__(self, spatial=False):
        super().__init__(Geometry(), spatial=spatial)

    def depth_of_a(self, a):
        a = np.maximum(np.asarray(a, float), 1e-12)
        return a*a/(2*self.g.radius)*(np.log(4*self.reference/a)+.5)

    def depth_derivative(self, a):
        a = np.maximum(np.asarray(a, float), 1e-12)
        return a/self.g.radius*np.log(4*self.reference/a)

    def a_of_depth(self, d):
        d = np.maximum(np.asarray(d, float), 0.)
        a = np.clip(np.sqrt(self.g.radius*d), 1e-8, .08)
        for _ in range(12):
            if np.all(np.abs(self.depth_of_a(a)-d)<=1e-16):
                break
            h = (self.depth_of_a(a)-d)/np.maximum(self.depth_derivative(a), 1e-5)
            a = np.clip(a-np.clip(h, -.5*a, .5*a), 1e-9, .08)
        return a

    def initial_log_scale(self, force, stiffness):
        a = .01
        for _ in range(15):
            a = math.sqrt(2*self.g.radius*force/stiffness/math.log(4*self.reference/a))
            a = float(np.clip(a, 1e-6, .05))
        return math.log(force/(a*a))

    def _at_a(self, x, a, theta):
        x, a = np.broadcast_arrays(np.asarray(x, float), np.asarray(a, float))
        c = np.exp(theta[0]); A = self.material(x[..., None]+a[..., None]*self.u, theta)
        f = c*a*a*(A@self.w)
        m = c*a**3*(A@(self.w*self.u))
        Ak = self.material(x[..., None]+a[..., None]*self.uk, theta)
        da = c*a*2*(Ak@self.wk)  # sum(wk)=1 at n=1/2
        return f, m, da

    def response(self, x, d, theta, stiffness=True):
        a = self.a_of_depth(d)
        f, m, da = self._at_a(x, a, theta)
        k = da/np.maximum(self.depth_derivative(a), 1e-12) if stiffness else None
        return np.where(np.asarray(d)>0,f,0.), np.where(np.asarray(d)>0,m,0.), k

    def inverse(self, x, f, theta):
        f = np.maximum(np.asarray(f, float), .005)
        a = np.clip(np.sqrt(f/np.exp(theta[0])), 1e-8, .08)
        for _ in range(16):
            ff, _, da = self._at_a(x, a, theta)
            if np.all(np.abs(ff-f)<=1e-13*np.maximum(f,1.)):
                break
            h = (ff-f)/np.maximum(da, 1e-6)
            a = np.clip(a-np.clip(h, -.6*a, .6*a), 1e-9, .08)
        return self.depth_of_a(a)
