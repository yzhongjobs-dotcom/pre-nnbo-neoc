"""Beber ICRA 2024 Fig. 2 from its vector data. Run after: pdftocairo -f 4 -l 4 -svg literature/papers/beber.pdf p4.svg, in the folder holding p4.svg."""
import re, json
import numpy as np
s = open('p4.svg', encoding='utf-8').read()
paths = re.findall(r'<path[^>]*?>', s)
def coords(p):
    d = re.search(r' d="([^"]+)"', p).group(1)
    nums = np.array(re.findall(r'-?\d+\.?\d*', d), float)
    return nums.reshape(-1, 2)
curves = {}
for p in paths:
    if ' d="' not in p or 'stroke:rgb' not in p:
        continue
    d = re.search(r' d="([^"]+)"', p).group(1)
    if d.count('L') < 1000:
        continue
    col = re.search(r'stroke:rgb\(([^)]+)\)', p).group(1)
    tr = re.search(r'matrix\(([^)]+)\)', p).group(1)
    if not tr.startswith('0.75'):      # the inset copies are rotated and magnified
        continue
    curves[col[:6]] = coords(p)
print({k: v.shape for k, v in curves.items()})
# axis frame: look for short paths with the same transform
frame = []
for p in paths:
    if ' d="' in p and 'stroke:rgb' in p and 'matrix(0.75' in p:
        d = re.search(r' d="([^"]+)"', p).group(1)
        if d.count('L') <= 5 and 'rgb(0%,0%,0%)' in p or d.count('L') <= 5 and 'rgb(14.' in p:
            c = coords(p)
            frame.append((d[:70], p[p.find('stroke:rgb'):p.find('stroke:rgb')+40]))
for f in frame[:30]:
    print(f)
b = curves[[k for k in curves if k.startswith('0%,44')][0]]
print('blue x range', b[:, 0].min(), b[:, 0].max(), 'y range', b[:, 1].min(), b[:, 1].max())
print('start', b[0], 'end', b[-1])
imax = int(np.argmax(b[:, 0])); print('x max index', imax, b[imax], 'y min index', int(np.argmin(b[:, 1])))
for i in np.linspace(0, len(b)-1, 25).astype(int):
    print(i, b[i].round(2))
np.save('blue_raw.npy', b)
import numpy as np, math
b = np.load('blue_raw.npy')
ex = (b[:, 0]-49.394674)/(340.998841-49.394674)*15.
fy = (121.455836-b[:, 1])/(121.455836-10.997503)*15.
n = len(ex)
print('points', n, 'peak penetration', ex.max().round(2), 'mm at', int(np.argmax(ex)), 'peak force', fy.max().round(2))
# smooth in index order
w = 41
sm = lambda x: np.convolve(np.pad(x, w//2, mode='edge'), np.ones(w)/w, mode='valid')
e, f = sm(ex), sm(fy)
top = int(np.argmax(e))
print('smoothed peak', e[top].round(2), f[top].round(2))
print('\nindex: penetration mm, force N')
for i in list(range(0, n, 400)):
    print(i, e[i].round(2), f[i].round(2))
# load and unload branches over penetration
load_end = int(np.argmax(e > e.max()-.3))
unload_start = top
le, lf = e[:load_end], f[:load_end]
ue, uf = e[unload_start:][::-1], f[unload_start:][::-1]
print('\nload ends at index', load_end, 'unload from', unload_start)
grid = np.arange(1., 14.1, 1.)
order = np.argsort(le); order2 = np.argsort(ue)
fl = np.interp(grid, le[order], lf[order]); fu = np.interp(grid, ue[order2], uf[order2])
print('e mm:', grid); print('load:', fl.round(2)); print('unload:', fu.round(2)); print('gap:', (fl-fu).round(2))
area = np.trapz(fl-fu, grid); wl = np.trapz(fl, grid)
print('loop area', round(area, 2), 'N mm, loading work', round(wl, 2), 'ratio', round(area/wl, 3), '-> equivalent sin(phi) about', round(area/wl/(2*math.pi)*2, 3), 'phase deg', round(math.degrees(math.asin(area/wl/math.pi)), 2))
# hold: force at top of the curve while penetration stays within .3 mm of its maximum
hold = np.where(e > e.max()-.3)[0]
print('\nhold samples', hold[0], hold[-1], 'force at start', f[hold[0]:hold[0]+100].mean().round(2), 'end', f[hold[-1]-100:hold[-1]].mean().round(2), 'max', f[hold].max().round(2), 'min', f[hold].min().round(2), 'penetration change', (e[hold[-1]]-e[hold[0]]).round(3))
# power-law fit of the load
m = (le > 1.) & (le < 14.)
p = np.polyfit(np.log(le[m]), np.log(lf[m]), 1); print('load exponent', p[0].round(3), 'coefficient F at 1 mm', np.exp(p[1]).round(3))
k = np.gradient(fl, grid); print('load tangent N/mm', k.round(2))
