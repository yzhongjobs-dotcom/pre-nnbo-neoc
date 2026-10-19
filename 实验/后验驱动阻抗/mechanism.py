"""Two-decision Bayesian impedance probe; no entropy reward or policy bank.

Known finite contact hypotheses, local prestressed linear dynamics and exact
encoders. This is a mechanism experiment, not B0 or a continuous scan benchmark.
"""
from itertools import product
from pathlib import Path
import hashlib
import json
import platform
import time

import numpy as np
import scipy
from scipy.linalg import expm, solve_continuous_lyapunov
from scipy.stats import qmc, norm
from numba import njit
from numpy.polynomial.hermite import hermgauss
from numpy.polynomial.legendre import leggauss

ROOT = Path(__file__).resolve().parent
SCALE = np.array([.001, .01])
MASS = .5
MODES = np.array(list(product([300., 1500.], [1.5, 18.])))
GAINS = np.array(list(product([0., 50., 150., 400., 1000.], [4., 20., 60.])))


def sample(n, seed):
    rng = np.random.default_rng(seed)
    return [dict(prior=rng.dirichlet(np.ones(4)*2).tolist(),
                 x0=[float(rng.normal(0, .05)), float(rng.normal(0, .05))],
                 speed=float(rng.uniform(.002, .008)),
                 probe=float(rng.uniform(.02, .06)),
                 future=float(rng.uniform(.3, .7)),
                 sigma=float(rng.choice([.015, .03, .06])),
                 force_noise=float(rng.choice([.03, .1, .3]))) for _ in range(n)]


def kernel(k, b, K, D, duration, r0, speed, sigma, order=32):
    """Physical cost as z.T S z+c; scaled x=[q/1mm,v/(1cm/s)]."""
    A = np.array([[0., 1.], [-(k+K)/MASS, -(b+D)/MASS]])
    G = np.array([[0., 0.], [0., (sigma/MASS)**2]])
    Cinf = solve_continuous_lyapunov(A, -G)
    vals, vecs = np.linalg.eig(A.astype(complex))
    inv = np.linalg.inv(vecs)
    def phi(t):
        return ((vecs*np.exp(vals*t))@inv).real
    a = K*speed/(k+K)
    c = (K*r0+D*speed-(b+D)*a)/(k+K)
    equilibrium = np.array([c, a])
    def transition(t):
        P = phi(t)
        beta = np.array([a*t+c, a])-P@equilibrium
        return P, beta
    P, beta = transition(duration)
    cov = Cinf-P@Cinf@P.T
    cov = (cov+cov.T)/2
    nodes, weights = leggauss(order)
    S = np.zeros((3, 3)); C = np.zeros((3, 3))
    for t, w in zip((nodes+1)*duration/2, weights*duration/2):
        E, shift = transition(t)
        T = np.eye(3); T[:2, :2] = E*SCALE[None, :]; T[:2, 2] = shift
        r = r0+speed*t
        force = np.array([k, b, 0.])/.3
        position = np.array([1., 0., -r])/.0015
        velocity = .2*np.array([0., 1., -speed])/.03
        command = np.sqrt(.015)*np.array([-K, -D, K*r+D*speed])/3.
        Q = np.outer(force, force)+np.outer(position, position)
        Q += np.outer(velocity, velocity)+np.outer(command, command)
        S += w*(T.T@Q@T)
        C += w*Q
    # Q's top-left 2x2 block is constant for this prescribed linear reference.
    Q2 = C[:2, :2]/duration
    O = solve_continuous_lyapunov(A.T, -Q2)
    S2 = O-P.T@O@P
    noise_cost = duration*np.trace(Q2@Cinf)-np.trace(S2@Cinf)
    scaled_cov = cov/np.outer(SCALE, SCALE)
    eigmin = np.linalg.eigvalsh(scaled_cov).min()
    assert eigmin > 0 and noise_cost >= -1e-10
    scaled_transition = P*SCALE[None, :]/SCALE[:, None]
    return S, max(0., noise_cost), scaled_transition, beta/SCALE, scaled_cov


def prepare(context, gains=GAINS, modes=MODES, order=32):
    nm, na = len(modes), len(gains)
    S = np.empty((2, nm, na, 3, 3)); c = np.empty((2, nm, na))
    P = np.empty((nm, na, 2, 2)); beta = np.empty((nm, na, 2))
    cov = np.empty((nm, na, 2, 2))
    for i, (k, b) in enumerate(modes):
        for a, (K, D) in enumerate(gains):
            for stage in range(2):
                h = context['probe'] if stage == 0 else context['future']
                r = 0. if stage == 0 else context['speed']*context['probe']
                out = kernel(k, b, K, D, h, r, context['speed'], context['sigma'], order)
                S[stage, i, a], c[stage, i, a] = out[:2]
                if stage == 0: P[i, a], beta[i, a], cov[i, a] = out[2:]
    return S, c, P, beta, cov


def quadrature(order):
    z, w = hermgauss(order); z *= np.sqrt(2); w /= np.sqrt(np.pi)
    indices = np.array(list(product(range(order), repeat=3)))
    return z[indices], np.prod(w[indices], axis=1)


def sobol_nodes(power=10, seed=11001):
    u = qmc.Sobol(3, scramble=True, seed=seed).random_base2(power)
    return norm.ppf(u), np.full(len(u), 1/len(u))


@njit(cache=True)
def costs(x, S, c):
    out = c.copy()
    for i in range(S.shape[0]):
        for a in range(S.shape[1]):
            H = S[i, a]
            out[i, a] += (H[0, 0]*x[0]**2+2*H[0, 1]*x[0]*x[1]
                         +H[1, 1]*x[1]**2+2*H[0, 2]*x[0]
                         +2*H[1, 2]*x[1]+H[2, 2])
    return out


@njit(cache=True)
def posterior(x, force, prior, means, inverse, logdet, force_noise, mode=2):
    nm = len(prior); logp = np.zeros(nm)
    for i in range(nm):
        dx = x-means[i]
        logp[i] = np.log(prior[i])-.5*(dx@inverse[i]@dx+logdet[i])
        if mode == 2:
            predicted = MODES[i, 0]*x[0]*SCALE[0]+MODES[i, 1]*x[1]*SCALE[1]
            logp[i] -= .5*((force-predicted)/force_noise)**2
    logp -= logp.max(); p = np.exp(logp); p /= p.sum()
    return p


def predictions(x, P, beta, cov):
    means = np.einsum('iajk,k->iaj', P, x)+beta
    return means, np.linalg.cholesky(cov), np.linalg.inv(cov), np.linalg.slogdet(cov)[1]


@njit(cache=True)
def scores(x0, prior, S, c, means, chol, inverse, logdet, force_noise, nodes, weights, mode):
    """mode0 freeze belief; mode1 future encoder update; mode2 full Bayes."""
    nm, na = means.shape[:2]; q = costs(x0, S[0], c[0])
    result = np.zeros(na)
    for a in range(na):
        result[a] = np.dot(prior, q[:, a])
        for i in range(nm):
            for j in range(len(weights)):
                x = means[i, a]+chol[i, a]@nodes[j, :2]
                force = MODES[i, 0]*x[0]*SCALE[0]+MODES[i, 1]*x[1]*SCALE[1]
                force += force_noise*nodes[j, 2]
                p = prior if mode == 0 else posterior(x, force, prior, means[:, a], inverse[:, a], logdet[:, a], force_noise, mode)
                value = costs(x, S[1], c[1])
                best = 1e100
                for b in range(na): best = min(best, np.dot(p, value[:, b]))
                result[a] += prior[i]*weights[j]*best
    return result


def open_loop_values(context, prepared):
    S, c, P, beta, cov = prepared
    x = np.asarray(context['x0']); prior = np.asarray(context['prior'])
    means = np.einsum('iajk,k->iaj', P, x)+beta
    q = costs(x, S[0], c[0]); out = np.zeros((len(GAINS), len(GAINS)))
    for i in range(len(MODES)):
        for a in range(len(GAINS)):
            y = np.r_[means[i, a], 1.]
            second = np.einsum('j,bjk,k->b', y, S[1, i], y)+c[1, i]
            second += np.einsum('bjk,kj->b', S[1, i, :, :2, :2], cov[i, a])
            out[a] += prior[i]*(q[i, a]+second)
    return out


@njit(cache=True)
def evaluate(x0, prior, S, c, means, chol, inverse, logdet, force_noise, samples, first_actions, fixed_second):
    q = costs(x0, S[0], c[0]); nmethod = len(first_actions)
    result = np.zeros(nmethod)
    # Conditional expectations of the full second block, MC over the FIRST
    # observed transition. Exact expected first-block cost is Rao-Blackwellized.
    for method in range(nmethod):
        a = first_actions[method]
        result[method] = np.dot(prior, q[:, a])
        for i in range(len(prior)):
            for j in range(len(samples)):
                x = means[i, a]+chol[i, a]@samples[j, :2]
                force = MODES[i, 0]*x[0]*SCALE[0]+MODES[i, 1]*x[1]*SCALE[1]
                force += force_noise*samples[j, 2]
                p = posterior(x, force, prior, means[:, a], inverse[:, a], logdet[:, a], force_noise)
                values = costs(x, S[1], c[1])
                b = fixed_second if method == 0 else np.argmin(p@values)
                result[method] += prior[i]*values[i, b]/len(samples)
    return result


def run_case(context, fixed, order=5, nmc=1024, seed=0, nodes_weights=None):
    tick = time.perf_counter(); prepared = prepare(context)
    setup_ms = (time.perf_counter()-tick)*1000
    S, c, P, beta, cov = prepared
    x = np.array(context['x0']); prior = np.array(context['prior'])
    means, chol, inv, logdet = predictions(x, P, beta, cov)
    nodes, weights = quadrature(order) if nodes_weights is None else nodes_weights
    actions = [fixed[0]]; values=[]; timing=[]
    for mode in range(3):
        tick = time.perf_counter()
        value = scores(x, prior, S, c, means, chol, inv, logdet, context['force_noise'], nodes, weights, mode)
        timing.append((time.perf_counter()-tick)*1000)
        actions.append(int(np.argmin(value))); values.append(value)
    samples = np.random.default_rng(seed).normal(size=(nmc, 3))
    result = evaluate(x, prior, S, c, means, chol, inv, logdet, context['force_noise'], samples, np.array(actions), fixed[1])
    return result, actions, values, timing, setup_ms


def pilot():
    cases = sample(8, 10401)
    prepared = [prepare(ctx) for ctx in cases]
    avg = np.mean([open_loop_values(ctx, pp) for ctx, pp in zip(cases, prepared)], axis=0)
    fixed = tuple(map(int, np.unravel_index(np.argmin(avg), avg.shape)))
    rows=[]
    for i, ctx in enumerate(cases):
        result, actions, values, times, setup = run_case(ctx, fixed, 5, 2048, 10402+i)
        rows.append(dict(J=result.tolist(), actions=actions, planning_ms=times, setup_ms=setup))
        print(i, result, actions, flush=True)
    output=dict(scope='Exploratory training/pilot only; not independent performance evidence',
                cases=cases, fixed=fixed, gains=GAINS.tolist(), modes=MODES.tolist(), rows=rows,
                source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                packages=dict(numpy=np.__version__,scipy=scipy.__version__),platform=platform.platform())
    (ROOT/'pilot.json').write_text(json.dumps(output, indent=2)+'\n')
    print('MEANS',np.mean([x['J'] for x in rows],axis=0),flush=True)


if __name__ == '__main__':
    pilot()
