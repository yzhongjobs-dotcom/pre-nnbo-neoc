import sys
sys.path.insert(0, r'E:\github_project\scan-impedance-update\experiments\rossi_failure_study')
from concurrent.futures import ProcessPoolExecutor
import scenes3 as S3, physics_plant3 as P3, method3 as M

O = dict(slide='long', pose_error=0., rho=0., force=12., dynamics='neural_white', phase=14., probe='convex', contact_law='foundation')
CH = {'as drawn': {}, 'offset 0': dict(sensor_delay_ms=0., delay_error=0.), 'offset exact': dict(delay_error=0.), 'mass exact': dict(mass_error=0.),
      'lag exact': dict(tau_error=0.), 'all exact': dict(mass_error=0., tau_error=0., delay_error=0.), 'offset 1.88, rest exact': dict(mass_error=0., tau_error=0., delay_error=0., sensor_delay_ms=1.88),
      'offset 1.0, rest as drawn': dict(sensor_delay_ms=1.)}


def job(n):
    sc = S3.cell_scene3('shift_fast', 'relaxed', 'nominal', 1, base=470900, overrides=dict(O, **CH[n]))
    pol = M.build('L', sc); log = P3.episode(sc, pol, 801); m = M.metric(log, pol, sc)
    return n, m['rmse'], m['accept'], m['est_abs'], sc['sensor_delay_ms'], sc['sensor_delay_hat'], sc['tool_mass_hat']/sc['tool_mass']-1, sc['tau_a_hat']/sc['tau_a_true']-1


if __name__ == '__main__':
    with ProcessPoolExecutor(8) as ex:
        for r in ex.map(job, CH):
            print('%-26s rmse %.3f accept %.2f gain err %.3f | delay %.2f hat %.2f mass err %+.3f lag err %+.3f' % r)
