"""Evidence-weighted bank of linear contact filters with different memories.

Observation, pivoted at the current encoder position (D is the force if the probe stayed there):
    z_j = k Y_j + W0 D_j + W1 r_j + W2 a_j + eps_j
    D+ = D + k dy + r dt + a dt^2/2,   r+ = r + a dt,   a+ = a
The motion prior is white jerk and must be broad: a prior narrower than the actual motion makes
slow tracking motion read as a soft contact. Each member differs only in how fast ln k may wander.
With slope=True the state also carries dk/dt, so a member predicts through a gain ramp.

Members never mix states. Weights come from each member's own predictive likelihood, discounted
over a horizon, so evidence accumulates even though one interval carries almost no information.
"""
import math
import numpy as np

K_MIN, K_MAX = 200., 3e5
TC = .008
A3 = np.array([[1., TC, TC**2/2], [0., 1., TC], [0., 0., 1.]])
Q3 = np.array([[TC**5/20, TC**4/8, TC**3/6], [TC**4/8, TC**3/3, TC**2/2], [TC**3/6, TC**2/2, TC]])


class Member:
    def __init__(self, q_lnk, q_motion, R, slope):
        self.ql, self.q, self.R, self.slope = q_lnk, q_motion, R, slope
        self.dim = 5 if slope else 4
        self.x = None; self.P = None; self.loglik = 0.
        self.iD = 2 if slope else 1

    def start(self, k, D, rel):
        d = self.iD
        self.x = np.zeros(self.dim); self.x[0] = k; self.x[d] = D
        P = np.zeros((self.dim, self.dim)); P[0, 0] = (rel*k)**2
        if self.slope:
            P[1, 1] = (k*self.ql)**2
        P[d, d] = 1.; P[d+1, d+1] = 1e2; P[d+2, d+2] = 1e4
        self.P = P; self.loglik = 0.

    def step(self, z, H4, dy):
        d = self.iD
        H = np.zeros(self.dim); H[0] = H4[0]; H[d:] = H4[1:]
        S = float(H@self.P@H+self.R)
        v = z-float(H@self.x)
        lp = -.5*(math.log(2*math.pi*S)+v*v/S)
        K = self.P@H/S
        x = self.x+K*v
        I = np.eye(self.dim)-np.outer(K, H)
        P = I@self.P@I.T+self.R*np.outer(K, K)
        F = np.eye(self.dim); F[d:, d:] = A3; F[d, 0] = dy
        kref = float(np.clip(x[0], K_MIN, K_MAX))
        Q = np.zeros((self.dim, self.dim)); Q[d:, d:] = kref*kref*self.q*Q3
        if self.slope:
            F[0, 1] = TC
            q = (kref*self.ql)**2
            Q[0, 0] = q*TC**3/3; Q[0, 1] = Q[1, 0] = q*TC**2/2; Q[1, 1] = q*TC
        else:
            Q[0, 0] = (kref*self.ql)**2*TC
        self.x = F@x; self.P = F@P@F.T+Q; self.P = (self.P+self.P.T)*.5
        self.x[0] = float(np.clip(self.x[0], K_MIN, K_MAX))
        return lp


class MemoryBank:
    def __init__(self, memories=(.1, .3, 1., 3.), q_motion=1.86, noise=.04, horizon=.25, slope=False, floor=1e-3):
        self.members = [Member(m, q_motion, noise**2, slope) for m in memories]
        self.lam = math.exp(-TC/horizon); self.floor = floor
        self.w = np.full(len(memories), 1./len(memories)); self.ready = False; self.lost = False

    def start(self, k, D, rel=.5):
        for m in self.members:
            m.start(k, D, rel)
        self.ready = True; self.lost = False

    def contact_lost(self):
        self.lost = True

    def step(self, z, H, dy):
        if self.lost:
            # Contact returned: force state restarts at the measured level, gain keeps its belief.
            for m in self.members:
                d = m.iD
                m.x[d] = z/max(H[1], 1e-6); m.x[d+1:] = 0.
                m.P[d:, :] = 0.; m.P[:, d:] = 0.
                m.P[d, d] = 1.; m.P[d+1, d+1] = 1e2; m.P[d+2, d+2] = 1e4
                m.loglik = 0.
            self.lost = False
        for m in self.members:
            m.loglik = self.lam*m.loglik+m.step(z, H, dy)
        ll = np.array([m.loglik for m in self.members])
        w = np.exp(ll-ll.max()); w = self.floor+(1-self.floor*len(w))*w/w.sum()
        self.w = w/w.sum()

    @property
    def k(self):
        return float(self.w@np.array([m.x[0] for m in self.members]))

    @property
    def sigma(self):
        ks = np.array([m.x[0] for m in self.members]); vs = np.array([m.P[0, 0] for m in self.members])
        mean = self.w@ks
        return math.sqrt(max(float(self.w@(vs+(ks-mean)**2)), 0.))
