"""Finite-footprint foundation / halfspace plant; no Method dependency.

SI units. Log columns 0:20 retain plant2 meanings; column 7 is E*, never q.
Column 5 is the sensor reading, including the force that accelerates the probe; column 18 is the true axial contact force.
Columns 20:26: Re(K*), |K*|, phase (rad), contact length, pressure centre,
minimum local cover (infinity off bone). K* includes robot series compliance,
but excludes support motion, actuator and sensor transfer functions.
"""
import math
import numpy as np
from numba import njit
from scipy.interpolate import PchipInterpolator

PLATFORM = dict(Tc=.008, ndel=1, tau_a=.019, fc=30., VMAX=.25)
LOGC, NS = 26, 37
GX, GW = np.polynomial.legendre.leggauss(4)
BX, BW = np.polynomial.legendre.leggauss(8)


@njit(cache=True)
def _lookup(field, dx, x):
    u = max(0., min(x/dx, len(field)-1.000000001))
    i = int(u); f = u-i
    return field[i]*(1.-f)+field[i+1]*f


@njit(cache=True)
def _bounds(d, theta, curvature, length):
    b = math.tan(theta); half = .5*length
    if curvature > 1.e-12:
        disc = b*b+4.*curvature*d
        if disc <= 0.:
            return 0., 0.
        root = math.sqrt(disc)
        return max(-half, (-b-root)/(2.*curvature)), min(half, (-b+root)/(2.*curvature))
    if abs(b) < 1.e-14:
        return (-half, half) if d > 0. else (0., 0.)
    return (-half, min(half, d/b)) if b > 0. else (max(-half, d/b), half)


@njit(cache=True)
def _rib_weight(x, ribs):
    w = 0.; cover = math.inf
    for row in ribs:
        center, width, edge, h = row
        distance = abs(x-center)
        a = (width*.5+edge*.5-distance)/edge
        if a > 0.:
            ww = 1. if a >= 1. else .5-.5*math.cos(math.pi*a)
            if ww > w:
                w = ww; cover = h
    return w, cover


@njit(cache=True)
def _hard(delta, v, invh, cap):
    """Pressure and tangent of a cover over bone. It hardens as the cover is used up until its stiffness per area reaches that of the rib and
    cage spread over the probe, then continues linearly, so the probe on bone is bounded by arm and rib in series and the law has no singularity."""
    if invh <= 0.:
        return v*delta, v
    uc = 1.-math.sqrt(v/cap) if cap > v else 0.
    u = delta*invh
    if u <= uc:
        return v*delta/(1.-u), v/(1.-u)**2
    dc = uc/invh
    return v*dc/(1.-uc)+cap*(delta-dc), cap


@njit(cache=True)
def _local(x, delta, c, qfield, efield, ribs):
    q = c[18]*_lookup(qfield, c[3], x-c[4]); e = c[19]*_lookup(efield, c[3], x-c[4])
    w, h = _rib_weight(x, ribs)
    if c[8] > 0.:
        w = 1.; h = c[8]
    v = q if c[5] == 0. else e
    p, stiffness = _hard(delta, v, w/h, c[13] if c[5] == 0. else c[14])
    return p, stiffness, h/w if w > 0. else math.inf


@njit(cache=True)
def _asym(coord, vc, span, c):
    if c[11] == 0. or vc == 0.:
        return 1.
    a0 = .5*span
    a1 = max(.2*a0, a0-c[9]*abs(vc))
    a2 = min(.5*span, a0+c[10]*abs(vc))
    return 1.+(a2-a1)/span*math.tanh(coord*math.copysign(1., vc)/.001)


@njit(cache=True)
def _section(xc, d, theta, vc, c, qfield, efield, ribs):
    length, width, radius = c[:3]
    curvature = .5/radius+.5/c[7]
    lo, hi = _bounds(d, theta, curvature, length)
    if hi <= lo:
        return 0., 0., 0., 0., math.inf, 0.
    area = width*(hi-lo)
    f = 0.; k = 0.; moment = 0.; pressure_moment = 0.; hmin = math.inf
    # Composite quadrature resolves 1 mm rib edges without a 2D dense mesh.
    thin = False
    for row in ribs:
        thin = thin or row[3] < .007
    nx = max(1, int(math.ceil((hi-lo)/(.001 if thin else .004))))
    xnodes = BX if thin else GX
    xweights = BW if thin else GW
    ny = 1 if c[6] == 0. else max(1, int(math.ceil(width/.004)))
    for sx in range(nx):
        mid = lo+(sx+.5)*(hi-lo)/nx
        hw = (hi-lo)/(2.*nx)
        for ix in range(len(xnodes)):
            xi = mid+hw*xnodes[ix]
            delta = max(d-curvature*xi*xi-xi*math.tan(theta), 0.)
            for sy in range(ny):
                # In long-axis scans the field is constant across the short axis.
                nq = 1 if c[6] == 0. else 4
                for iy in range(nq):
                    eta = 0. if nq == 1 else -.5*width+(sy+.5)*width/ny+.5*width/ny*GX[iy]
                    wt = hw*xweights[ix]*(width if nq == 1 else .5*width/ny*GW[iy])
                    asym = _asym(xi if c[6] == 0. else eta, vc, length if c[6] == 0. else width, c)
                    pos = xc+(xi if c[6] == 0. else eta)
                    if c[12] > 0.:
                        p = 0.; stiffness = 0.; h = math.inf
                        for transverse in range(int(c[17])):
                            jj = 26+4*transverse
                            invh = c[jj+2]
                            if invh > 0.:
                                h = min(h, 1./invh)
                            pi, si = _hard(delta, c[jj if c[5] == 0. else jj+1], invh, c[jj+3]*(c[13] if c[5] == 0. else c[14]))
                            p += pi; stiffness += si
                    else:
                        p, stiffness, h = _local(pos, delta, c, qfield, efield, ribs)
                    hmin = min(hmin, h)
                    f += wt*p*asym
                    pressure_moment += wt*p*xi*asym
                    k += wt*stiffness*asym
                    moment += wt*(p if c[5] == 0. else stiffness)*xi*asym
    if c[5] != 0.:
        # Same local tangent-hardening factor as the foundation; force is its integral in d.
        scale = 1.1/math.sqrt(area)
        k *= scale; moment *= scale
    return f, k, moment, hi-lo, hmin, pressure_moment


@njit(cache=True)
def _contact_integral(xc, d, theta, vc, c, qfield, efield, ribs):
    f, k, moment, length, h, pressure_moment = _section(xc, d, theta, vc, c, qfield, efield, ribs)
    if f <= 0. or length == 0.:
        return 0., k, 0., length, h
    if c[5] == 0. or f >= 1.e20:
        return f, k, moment/f, length, h
    # Halfspace: the force whose tangent is k = 1.1 E sqrt(A). While the contact grows a convex face has F = 1.2 (1.1/sqrt(A)) F_E
    # with F_E the same integral on E, exactly when aligned and homogeneous; once full, it continues at constant k.
    full = (.5/c[2]+.5/c[7])*(.5*c[0])**2+abs(math.tan(theta))*.5*c[0]
    scale = 1.1/math.sqrt(c[1]*length)
    if d <= full:
        return 1.2*scale*f, k, pressure_moment/f, length, h
    cap = 0.
    if full > 0.:
        cap = _section(xc, full, theta, vc, c, qfield, efield, ribs)[0]
    return scale*(f+.2*cap), k, pressure_moment/f, length, h


@njit(cache=True)
def _measured(f, k, c, active=False):
    """Monotone map from the reference force to the force the measured indentation curve gives on this probe; pose, ribs and field act through f."""
    n = int(c[16])
    if n == 0 or f >= 1.e20:
        return f, k
    i0 = int(c[24] if active else c[15])
    ref = c[i0:i0+n]; val = c[i0+n:i0+2*n]; slope = c[i0+2*n:i0+3*n]
    if f <= ref[0]:
        return val[0]+(f-ref[0])*slope[0], k*slope[0]
    if f >= ref[n-1]:
        return val[n-1]+slope[n-1]*(f-ref[n-1]), k*slope[n-1]
    j = int(np.searchsorted(ref, f))-1
    w = (f-ref[j])/(ref[j+1]-ref[j])
    h = ref[j+1]-ref[j]
    a = 2.*val[j]-2.*val[j+1]+h*(slope[j]+slope[j+1])
    b = -3.*val[j]+3.*val[j+1]-h*(2.*slope[j]+slope[j+1])
    value = ((a*w+b)*w+h*slope[j])*w+val[j]
    tangent = (3.*a*w*w+2.*b*w)/h+slope[j]
    return value, k*tangent


@njit(cache=True)
def _contact(xc, d, theta, vc, c, qfield, efield, ribs, activation=0.):
    f, k, cop, length, h = _reference(xc, d, theta, vc, c, qfield, efield, ribs)
    # the muscle moves the force between the same person's relaxed and activated curves; the activation time course is an assumption
    fa, ka = _measured(f, k, c, True)
    f, k = _measured(f, k, c)
    return c[20]*((1.-activation)*f+activation*fa), c[20]*((1.-activation)*k+activation*ka), cop, length, h


@njit(cache=True)
def _reference(xc, d, theta, vc, c, qfield, efield, ribs):
    if c[6] == 0.:
        return _contact_integral(xc, d, theta, vc, c, qfield, efield, ribs)
    # Transverse field samples are reused at every depth quadrature point.
    short = c.copy(); short[6] = 0.; short[11] = 0.; short[12] = 1.
    thin = False
    for row in ribs:
        thin = thin or row[3] < .007
    width = c[1]; n = max(1, int(math.ceil(width/(.001 if thin else .004))))
    num = 0
    for j in range(n):
        for k in range(4):
            eta = -.5*width+(j+.5)*width/n+.5*width/n*GX[k]
            weight = .5/n*GW[k]*_asym(eta, vc, width, c)
            x = xc+eta; w, h = _rib_weight(x, ribs)
            if c[8] > 0.:
                w = 1.; h = c[8]
            q = _lookup(qfield, c[3], x-c[4]); e = _lookup(efield, c[3], x-c[4])
            # Equal cover depths have identical rational pressure factors.
            invh = w/h; group = num
            for previous in range(num):
                if short[26+4*previous+2] == invh:
                    group = previous
                    break
            jj = 26+4*group
            short[jj] += weight*q*c[18]; short[jj+1] += weight*e*c[19]; short[jj+2] = invh; short[jj+3] += weight
            if group == num:
                num += 1
    short[17] = num
    return _contact_integral(xc, d, theta, vc, short, qfield, efield, ribs)


def contact(sc, d, x=None, theta=0., velocity=0., activation=0.):
    args = (float(sc['xcmd'][0] if x is None else x), float(d), float(theta), float(velocity), _geometry(sc), sc['q_field'], sc['E_field'], sc['ribs'], float(activation))
    f, k, cop, length, h = _contact(*args)
    if sc['contact_law'] == 'halfspace':
        # on a convex face over a cover the closed-form force only approximates the integral of k; its slope is the tangent
        eps = 1.e-7
        k = (_contact(args[0], args[1]+eps, *args[2:])[0]-_contact(args[0], max(args[1]-eps, 0.), *args[2:])[0])/(args[1]+eps-max(args[1]-eps, 0.))
    return f, k, cop, length, h


def _geometry(sc):
    if sc['contact_law'] not in ('foundation', 'halfspace'):
        raise ValueError(sc['contact_law'])
    if sc['slide'] not in ('long', 'short'):
        raise ValueError(sc['slide'])
    g = sc['geometry']
    base = np.array([g['length'], g['width'], g['radius'], sc['dx'], sc['field_x0'],
                     float(sc['contact_law'] == 'halfspace'), float(sc['slide'] == 'short'),
                     sc['surface_radius'], sc['cover'] if sc.get('full_cover') else 0.,
                     sc['k1'], sc['k2'], float(sc['asymmetry'] and sc['mu_k'] > 0.), 0.,
                     sc.get('k_rib', 1.e12)/(g['length']*g['width']), sc.get('k_rib', 1.e12)/(1.1*math.sqrt(g['length']*g['width'])), 0., 0., 0.,
                     sc['q_base'], sc['E_base'], sc.get('q_scale', 1.), 0., sc.get('k_rib', 1.e12), 0., 0., 0.])
    out = np.r_[base, np.zeros(16*int(math.ceil(g['width']/.001)))]
    table = sc.get('curve_table')
    if table is not None:
        active = sc.get('active_table', table)
        # the rib limit applies to the force after the map, at the steepest slope of either curve
        steepest = max(float(np.max(table[2])), float(np.max(active[2])))*sc.get('q_scale', 1.)
        out[13] /= steepest; out[14] /= steepest
        out[15] = len(out); out[16] = len(table[0])
        out = np.r_[out, table[0], table[1], table[2]]
        out[24] = len(out)
        out = np.r_[out, active[0], active[1], active[2]]
    if min(out[13], out[14], out[20], out[22]) <= 0.:
        raise ValueError('Nonpositive contact scale or cap')
    return out


REMUS_R = .010        # hemisphere radius of the measured curves, m


def _hemisphere_area(delta):
    d = np.minimum(delta, REMUS_R)
    return np.pi*(2.*REMUS_R*d-d*d)


def positive_domain(poly, fit_max_mm):
    """The rising stretch of a fitted curve that holds the 5 mm reference; outside it the curve continues along its end tangent."""
    roots = np.roots(np.polyder(poly))
    breaks = [0.]+sorted(float(r.real) for r in roots if abs(r.imag) < 1.e-8 and 0. < r.real < fit_max_mm)+[fit_max_mm]
    intervals = [(a, b) for a, b in zip(breaks[:-1], breaks[1:]) if np.polyval(np.polyder(poly), .5*(a+b)) > 0.]
    valid = [(a, b) for a, b in intervals if a < min(5., fit_max_mm) <= b]
    if len(valid) != 1:
        raise ValueError('No positive tangent interval containing the regularization point')
    a, b = valid[0]
    return a+1.e-6 if a > 0. else 0., b-1.e-6 if b < fit_max_mm else b


def probe_curve(poly, geometry, gamma, depths, fit_max_mm=None):
    """Force on the probe from a measured force-depth polynomial on the r = 10 mm hemisphere (Remus Table 4, N and mm, descending powers).

    Equal strain s = depth/sqrt(area) means equal material state; at equal s the tangent scales with contact area as area^gamma, .5 for a
    half-space, where stiffness follows sqrt(A), and 1 for a thin layer, where it follows A. Geometry alone therefore sets the shape on the probe.
    """
    L, W, R = geometry['length'], geometry['width'], geometry['radius']
    a = np.full_like(depths, L/2.) if math.isinf(R) else np.minimum(np.sqrt(2.*R*depths), L/2.)
    area = W*2.*a
    s = depths/np.sqrt(area)
    dh = np.geomspace(1.e-5, .05, 4000)
    sh = dh/np.sqrt(_hemisphere_area(dh))
    # The 5 mm material-state floor is a regularization assumption, not a measured threshold.
    delta = np.maximum(np.interp(s, sh, dh), .005)
    if fit_max_mm is None:
        raise ValueError('fitMaxDispl is required for measured curves')
    lower, upper = positive_domain(poly, fit_max_mm)
    k_h = np.polyval(np.polyder(poly), np.clip(delta*1000., lower, upper))*1000.
    k = k_h*(area/_hemisphere_area(delta))**gamma
    if not np.all(np.isfinite(k)) or np.any(k <= 0.):
        raise ValueError('Nonpositive transferred tangent')
    force = np.r_[0., np.cumsum(.5*(k[1:]+k[:-1])*np.diff(depths))]
    return force, k


def attach_curve(sc):
    """Table from the plant's reference force to the measured force for this scene's probe and transfer assumption."""
    g = sc['geometry']
    ref = dict(sc, ribs=np.empty((0, 4)), full_cover=False, surface_radius=math.inf, asymmetry=False, q_field=np.ones(2),
               E_field=np.ones(2), dx=1., field_x0=0., curve_table=None, active_table=None, q_scale=1.)
    depths = np.geomspace(1.e-6, .10, 400)
    fr = np.array([contact(ref, d, 0.)[0] for d in depths])
    fine = np.geomspace(1.e-8, .10, 4000)
    for name, poly, limit in (('curve_table', sc['remus_poly'], sc['fit_max_mm']), ('active_table', sc['active_poly'], sc['active_fit_max_mm'])):
        force, k = probe_curve(poly, g, sc['transfer_gamma'], fine, limit)
        force = force+k[0]*fine[0]/(1.5 if math.isfinite(g['radius']) else 1.)
        val = np.r_[0., np.interp(depths, fine, force)]
        x = np.r_[0., fr]
        slopes = PchipInterpolator(x, val).derivative()(x)
        if np.any(np.diff(x) <= 0.) or np.any(np.diff(val) <= 0.) or np.any(slopes <= 0.):
            raise ValueError('Nonmonotone measured map or nonpositive slope')
        sc[name] = (x, val, slopes)


def solve_depth(sc, force, x=None, theta=0.):
    if force <= 0.:
        raise ValueError('equilibrium force must be positive')
    geom = _geometry(sc)
    xx = float(sc['xcmd'][0] if x is None else x)
    curvature = .5/geom[2]+.5/geom[7]; b = math.tan(theta); half = .5*geom[0]
    vertex = min(half, max(-half, -b/(2.*curvature))) if curvature > 1.e-12 else (-half if b >= 0. else half)
    lo = curvature*vertex*vertex+b*vertex-1.e-12
    hi = max(.001, lo+.001)
    def value(depth):
        return _contact(xx, depth, theta, 0., geom, sc['q_field'], sc['E_field'], sc['ribs'])
    while value(hi)[0] < force:
        hi = lo+2.*(hi-lo)
        if hi > 1.:
            raise ValueError('no equilibrium within one metre')
    d = .5*(lo+hi)
    for iteration in range(200):
        f, k, _, _, _ = value(d)
        if abs(f-force) < 1.e-9 or hi-lo < 1.e-12:
            return d
        if f > force:
            hi = d
        else:
            lo = d
        # Near a thin cover the force is steep on one side, so Newton alone creeps; every third step halves the bracket.
        candidate = d-(f-force)/k if k > 0. and iteration % 3 != 2 else math.inf
        d = candidate if lo < candidate < hi else .5*(lo+hi)
    raise RuntimeError('equilibrium solve did not converge')


def compliance(sc):
    off = sc['rho']/math.sqrt(sc['kr']*sc['K_t'])
    return np.array([[1./sc['kr'], off], [off, 1./sc['K_t']]])


def calibration(sc):
    # A fixed, zero-slide normal reference on the uniform soft tissue: independent of FD, position, ribs and scan regime.
    ref = dict(sc, ribs=np.empty((0, 4)), full_cover=False, surface_radius=math.inf, asymmetry=False, q_field=np.ones(2),
               E_field=np.ones(2), dx=1., field_x0=0., q_scale=1.)
    d = solve_depth(ref, 8., 0.)
    k = contact(ref, d, 0.)[1]
    return math.tan(math.radians(sc['phi_hc']))*k/(2.*math.pi*20.8*8.)


def complex_gain(sc, k, force, frequency=20.8, theta=0., lam=None):
    w = 2.*math.pi*frequency
    memory = 1.+np.sum(sc['alpha']*(1j*w*sc['tau'])/(1.+1j*w*sc['tau']))
    kc = k*memory+1j*w*(calibration(sc) if lam is None else lam)*force
    return kc/(1.+kc*math.cos(theta)**2/sc['kr'])*math.cos(theta)**2


def hemisphere(sc, d, radius=.01):
    """Transfer to the finite Remus hemisphere, including travel beyond r.

    The halfspace circular-punch factor is 2 E* a, rather than the rectangular
    beta approximation. Hertz d<<r is a separate approximation to this shape.
    """
    d = max(float(d), 0.)
    r = radius
    if sc['contact_law'] == 'foundation':
        volume = math.pi*(r*d*d-d**3/3.) if d <= r else 2.*math.pi*r**3/3.+math.pi*r*r*(d-r)
        return sc['q_base']*volume
    # Integral of circular flat-punch stiffness 2 E a(s), a(s)=sqrt(2rs-s^2).
    s = min(d, r); z = s-r
    integral = .5*(z*math.sqrt(max(r*r-z*z, 0.))+r*r*math.asin(z/r))+.25*math.pi*r*r
    return 2.*sc['E_base']*(integral+r*max(d-r, 0.))


@njit(cache=True)
def _friction(fn, vc, length, c, geom):
    mu = c[13]+(c[14]-c[13])*math.exp(-(vc/c[15])**2)
    span = length if geom[6] == 0. else geom[1]
    a0 = .5*span
    a1 = max(.2*a0, a0-geom[9]*abs(vc)); a2 = min(.5*(geom[0] if geom[6] == 0. else geom[1]), a0+geom[10]*abs(vc))
    wedge = 0.
    if geom[11] > 0. and math.isfinite(geom[2]) and geom[6] == 0.:
        wedge = 4.*fn*(a2-a1)/(3.*math.pi*geom[2])*math.tanh(vc/c[16])
    return mu*fn*math.tanh(vc/c[16])+c[17]*vc+wedge if fn > 0. else 0.


@njit(cache=True)
def _theta(x, c):
    # Positive where the surface falls ahead of the motion: the pressure moves behind and friction adds to the axial force.
    return c[21]+(1.-c[22])*math.atan((x-c[28])/c[27])


@njit(cache=True)
def _tick(s, buf, cmd, c, geom, qfield, efield, ribs, hs, xcmd, u, fext, act, rip,
          alpha, tau, noise, kk, log, nlog):
    dt = .0005; Tc = c[31]; ndel = len(buf)-1
    cmd = min(max(cmd, buf[-1]-c[12]*Tc), buf[-1]+c[12]*Tc)
    for j in range(ndel):
        buf[j] = buf[j+1]
    buf[-1] = cmd
    ei = np.exp(-dt/tau); weights = alpha*tau/dt*(1.-ei)
    for j in range(int(round(Tc/dt))):
        t = s[0]; it = min(int(round(t/dt))//2, len(u)-2); fr = (t-it*.001)/.001
        uu = u[it]*(1.-fr)+u[it+1]*fr; ud = (u[it+1]-u[it])*1000.
        activation = min(1., max(0., act[it]*(1.-fr)+act[it+1]*fr))
        oldxs = s[1]
        s[1] += (1.-math.exp(-dt/c[5]))*(xcmd[it]*(1.-fr)+xcmd[it+1]*fr-s[1])
        vs = (s[1]-oldxs)/dt
        s[4] += dt*(c[2]*(uu-s[3])+c[1]*(ud-s[4])+fext[it]-s[16])/c[0]
        s[3] += dt*s[4]
        xold = s[7]; vold = s[8]; vc = vold; fn = s[14]; ftax = s[19]
        d = s[13]; fe = s[9]; ke = 0.; cop = 0.; length = 0.; cover = math.inf
        for outer in range(3):
            damping = 2.*c[25]*math.sqrt(c[23]*c[24])
            oldct = math.cos(s[17]); oldst = math.sin(s[17])
            ft0 = (_friction(fn, vc, length, c, geom)-c[17]*vc)*oldct-fn*oldst
            vc = (vold+dt/c[24]*(c[23]*(s[1]-xold-c[26]*s[16])+damping*vs-ft0))/(1.+dt/c[24]*(damping+c[17]*oldct+dt*c[23]))
            xc = xold+dt*vc; theta = _theta(xc, c)
            ct = math.cos(theta); st = math.sin(theta)
            h = _lookup(hs, geom[3], xc-geom[4])
            base = s[2]+s[3]+h+rip[it]
            cnn = 1./c[3]
            # The residual is monotone in d. Past the local cover the force is unbounded, so keep a bracket: the compliance in
            # series with the bone then sets the depth, and a probe sliding onto a thin rib stays well posed.
            below = -math.inf; above = math.inf; back = .0002
            for iteration in range(40):
                fe, ke, cop, length, cover = _contact(xc, d, theta, vc, geom, qfield, efield, ribs, activation)
                if fe >= 1.e20:
                    above = d
                    d = .5*(d+below) if math.isfinite(below) else d-back
                    back *= 2.
                    continue
                qsum = np.sum(s[10:13]*ei+weights*(fe-s[9]))
                rate = 1.+c[4]*(3.*d-4.*s[13]+s[22])/(2.*dt)
                fn = max((fe+qsum)*rate, 0.)
                ft = _friction(fn, vc, length, c, geom)
                fax = fn*ct+ft*st; ftax = ft*ct-fn*st
                # the arm carries the contact force and the force accelerating the probe
                residual = d-(base-(fax+c[33]*s[24])*cnn-c[26]*ftax)*ct
                if residual > 0.:
                    above = d
                else:
                    below = d
                derivative = ke*(1.+np.sum(weights))*max(rate, 0.)+1.5*(fe+qsum)*c[4]/dt
                slope = 1.+max(derivative, 0.)*ct*ct*cnn
                change = min(max(residual/slope, -.002), .002)
                if abs(change) < 1.e-12:
                    break
                d -= change
                if d >= above or d <= below:
                    if math.isfinite(above) and math.isfinite(below):
                        d = .5*(above+below)
                    else:
                        d = above-back if math.isfinite(above) else below+back
                        back *= 2.
        # Re-evaluate at the accepted implicit state before committing Prony history.
        fe, ke, cop, length, cover = _contact(xc, d, theta, vc, geom, qfield, efield, ribs, activation)
        if fe >= 1.e20 and math.isfinite(below):
            d = below
            fe, ke, cop, length, cover = _contact(xc, d, theta, vc, geom, qfield, efield, ribs, activation)
        s[10:13] = s[10:13]*ei+weights*(fe-s[9])
        normal_velocity = (3.*d-4.*s[13]+s[22])/(2.*dt)
        hc_factor = max(1.+c[4]*normal_velocity, 0.)
        hereditary_force = fe+np.sum(s[10:13])
        fn = max(hereditary_force*hc_factor, 0.)
        if length == 0.:
            s[10:13] = 0.; fn = 0.
        ft = _friction(fn, vc, length, c, geom)
        fax = fn*ct+ft*st; ftax = ft*ct-fn*st
        s[22] = s[13]
        s[7] = xc; s[8] = vc; s[9] = fe; s[13] = d; s[14] = fn; s[15] = ft; s[16] = fax; s[17] = theta
        # State 19 is reserved for sensor Ft below; do not use it as a true load.
        if int(round(t/dt)) % 2 == 0:
            k = min(kk, len(noise)-1); ns = c[7]
            s[6] += 2.e-4*ns*noise[k, 0]; s[18] += 2.e-4*ns*noise[k, 2]
            # The sensor sits between arm and probe, so it reads the contact force plus the force accelerating the probe (Beber eq. 3-4).
            # The servo model has a velocity jump at each command step, which no real probe follows, so acceleration is taken over one control period.
            velocity = (buf[0]-s[2])/c[5]
            s[24] = (velocity-s[29])/(Tc if Tc > 0. else .008)
            s[29:36] = s[30:37].copy(); s[36] = velocity
            s[28] = s[27]; s[27] = s[26]; s[26] = s[25]; s[25] = fax+c[33]*s[24]
            # time offset of the force stream against the encoder, fractions of a sample by interpolation
            lag = int(c[34]); frac = c[34]-lag
            sensed = (1.-frac)*s[25+lag]+frac*s[25+min(lag+1, 3)]
            s[20] += 2.e-4*(7.5e-4/.08)*ns*noise[k, 4]
            s[5] += c[6]*(sensed+s[6]+c[30]*ns*noise[k, 1]-s[5])
            s[19] += c[6]*(ftax+s[18]+c[30]*ns*noise[k, 3]-s[19])
            s[21] += c[6]*(fn*cop+ft*c[29]+s[20]+7.5e-4*(c[30]/.08)*ns*noise[k, 5]-s[21])
            kk += 1
        if j % int(round(c[32]/dt)) == 0:
            w = 2.*math.pi*20.8
            if geom[5] != 0. and fn > 0.:
                # On a convex face over a cover the closed-form force only approximates the integral of k; truth is its own slope.
                ke = (_contact(xc, d+1.e-6, theta, vc, geom, qfield, efield, ribs, activation)[0]-_contact(xc, d-1.e-6, theta, vc, geom, qfield, efield, ribs, activation)[0])/2.e-6
            tangent = ke*hc_factor if fn > 0. else 0.
            kc = complex(tangent, w*c[4]*max(hereditary_force, 0.)) if fn > 0. else 0.j
            kl = tangent; ki = tangent*(1.+np.sum(alpha))
            for q in range(3):
                z = complex(0., w*tau[q]); kc += tangent*alpha[q]*z/(1.+z)
                z = complex(0., 2.*math.pi*tau[q]); kl += tangent*alpha[q]*(z/(1.+z)).real
            friction_slope = _friction(1., vc, length, c, geom)-c[17]*vc
            axial = ct+friction_slope*st; lateral = friction_slope*ct-st
            series = ct*(axial*cnn+c[26]*lateral)
            kt = kc/(1.+kc*series)*axial*ct
            log[nlog, 0] = t; log[nlog, 1] = s[1]; log[nlog, 2] = s[2]; log[nlog, 3] = fn
            log[nlog, 4] = s[3]; log[nlog, 5] = s[5]; log[nlog, 6] = d
            log[nlog, 7] = geom[19]*_lookup(efield, geom[3], xc-geom[4]); log[nlog, 8] = kt.real
            log[nlog, 9] = kl/(1.+kl*series)*axial*ct; log[nlog, 10] = ki/(1.+ki*series)*axial*ct
            log[nlog, 11] = buf[-1]; log[nlog, 12] = s[19]; log[nlog, 13] = s[21]
            log[nlog, 14] = uu; log[nlog, 15] = act[it]; log[nlog, 16] = vc; log[nlog, 17] = theta
            log[nlog, 18] = fax; log[nlog, 19] = ft
            log[nlog, 20] = kt.real; log[nlog, 21] = abs(kt); log[nlog, 22] = math.atan2(kt.imag, kt.real)
            log[nlog, 23] = length; log[nlog, 24] = cop; log[nlog, 25] = cover
            nlog += 1
        s[2] += (1.-math.exp(-dt/c[5]))*(buf[0]-s[2]); s[0] = t+dt
        if not math.isfinite(fn) or fn > 5.e3 or abs(s[2]) > .5:
            return -nlog, kk
    return nlog, kk


def _constants(sc, pl):
    C = compliance(sc)
    lp = pl['fc'] if sc.get('sensor_lp') in (None, False) else (100. if sc['sensor_lp'] is True else float(sc['sensor_lp']))
    return np.array([sc['m'], sc['c'], sc['k'], sc['kr'], calibration(sc), pl['tau_a'],
                     1.-math.exp(-2.*math.pi*lp*.001), sc['noise_scale'], pl['ndel'], 0., 0., 0., pl['VMAX'],
                     sc['mu_k'], sc['mu_s'], sc['v_s'], sc['v_e'], sc['k_d'], 0., 0., 0., sc['theta0'],
                     float(sc['tracking']), sc['K_t'], sc['m_t'], sc['zeta_t'], C[0, 1], sc['surface_radius'],
                     sc['surface_center'], sc['l_tool'], sc['sensor_sigma'], pl['Tc'], pl.get('log_dt', pl['Tc']), sc.get('tool_mass', 0.), float(sc.get('sensor_delay_ms', 0))])


def _init(sc, pl):
    xs = float(sc['xcmd'][0]); x = xs; C = compliance(sc)
    for _ in range(12):
        theta = sc['theta0']+(not sc['tracking'])*math.atan((x-sc['surface_center'])/sc['surface_radius'])
        d = solve_depth(sc, sc['FD'], x, theta); f, ke, cop, _, _ = contact(sc, d, x, theta)
        fax = f*math.cos(theta); ftax = -f*math.sin(theta)
        target = xs-C[0, 1]*fax-C[1, 1]*ftax
        if abs(target-x) < 1.e-12:
            break
        x = target
    z = sc['u'][0]+(sc['fext'][0]-fax)/sc['k']
    h = np.interp(x, sc['field_x0']+np.arange(len(sc['h_s']))*sc['dx'], sc['h_s'])
    y = d/math.cos(theta)-z-h+C[0, 0]*fax+C[0, 1]*ftax
    s = np.zeros(NS); s[1] = xs; s[2] = y; s[3] = z; s[5] = fax
    s[7] = x; s[9] = f; s[13] = d; s[14] = f
    s[16] = fax; s[17] = theta; s[19] = ftax; s[21] = f*cop; s[22] = d; s[25:29] = fax
    gap = float(sc.get('start_gap', 0.))
    if gap > 0.:
        # starts free of the tissue, above the surface by gap
        y = float(sc['u'][0]+h*0.-gap-sc['h_s_at_start'] if 'h_s_at_start' in sc else sc['u'][0]-h-gap)
        s[2] = y; s[3] = sc['u'][0]; s[5] = 0.; s[9] = 0.; s[13] = -gap; s[14] = 0.; s[16] = 0.; s[19] = 0.; s[21] = 0.; s[22] = -gap; s[25:29] = 0.
    return s, np.full(int(pl['ndel'])+1, y)


class Obs:
    """Measured channels only; truth columns are masked even in .log."""
    def __init__(self, log, n):
        self.n = n
        self.log = np.full_like(log[:n], np.nan)
        cols = [0, 1, 2, 5, 11, 12, 13]
        self.log[:, cols] = log[:n, cols]

    t = property(lambda self: self.log[:, 0])
    x = property(lambda self: self.log[:, 1])
    y = property(lambda self: self.log[:, 2])
    force = property(lambda self: self.log[:, 5])
    y_cmd = property(lambda self: self.log[:, 11])
    ft = property(lambda self: self.log[:, 12])
    my = property(lambda self: self.log[:, 13])


def _noise(rng, n, sensor):
    if not sensor:
        return np.zeros((n, 6))
    return np.c_[rng.standard_normal((n, 2)), rng.standard_normal((n, 4))]


def _run(sc, seed, sensor, platform, commands=None, policy=None, diagnostic=False):
    pl = dict(PLATFORM if platform is None else platform)
    # The plant runs on the true system; the policy is told only what calibration can give it.
    calibrated = dict(pl)
    if 'tau_a_true' in sc and platform is None:
        pl['tau_a'] = sc['tau_a_true']; calibrated['tau_a'] = sc['tau_a_hat']
    if abs(pl['Tc']/.0005-round(pl['Tc']/.0005)) > 1.e-8:
        raise ValueError('Tc must be an integer multiple of 0.5 ms')
    nt = int(round(sc['T']/pl['Tc'])) if commands is None else len(commands)
    if nt*pl['Tc'] > sc['T']+1.e-9:
        raise ValueError('commands exceed scene drive duration')
    log_dt = pl.get('log_dt', pl['Tc'])
    if log_dt < .0005 or abs(log_dt/.0005-round(log_dt/.0005)) > 1.e-8 or abs(pl['Tc']/log_dt-round(pl['Tc']/log_dt)) > 1.e-8:
        raise ValueError('log_dt must divide Tc and be a multiple of 0.5 ms')
    log = np.zeros((int(round(nt*pl['Tc']/log_dt))+2, LOGC)); s, buf = _init(sc, pl); c = _constants(sc, pl)
    geom = _geometry(sc); noise = _noise(np.random.default_rng(seed), int(sc['T']*1000.)+10, sensor)
    if policy is not None:
        policy.reset(calibrated)
    n = 0; kk = 0
    for j in range(nt):
        seen = s[5]
        if policy is not None:
            log[n, [0, 1, 2, 5, 11, 12, 13]] = [s[0], s[1], s[2], s[5], buf[-1], s[19], s[21]]
            obs = Obs(log, n+1)
            if diagnostic and n > 0:
                obs.truth_k = float(log[n-1, 21]); obs.truth_re = float(log[n-1, 20])
            command = float(policy(obs))
        else:
            command = float(commands[j])
        if not np.isfinite(command):
            command = float(buf[-1])
        n0 = n
        n, kk = _tick(s, buf, command, c, geom, sc['q_field'], sc['E_field'], sc['ribs'], sc['h_s'],
                      sc['xcmd'], sc['u'], sc['fext'], sc['act'], sc['rip'], sc['alpha'], sc['tau'],
                      noise, kk, log, n)
        if n < 0:
            return None
        log[n0, 5] = seen
    return log[:n]


def replay(sc, commands, seed=0, sensor=False, platform=None):
    return _run(sc, seed, sensor, platform, commands=np.asarray(commands, float))


def episode(sc, policy, seed, platform=None, sensor=True, diagnostic=False):
    return _run(sc, seed, sensor, platform, policy=policy, diagnostic=diagnostic)
