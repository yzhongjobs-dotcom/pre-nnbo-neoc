"""Plot the motion-isolation boundary with the established Python backend.

Figure contract: do the frozen controllers maintain force while following a
moving homogeneous surface, and does the answer change on stiff material?
Quantitative grid: columns are the first soft/stiff push scene (indices 0/2),
rows show force and actual normal displacement. Four illustrative controllers;
all thirteen methods and all forty scenes remain in the numerical tables.
Neural model 0 is shown; no smoothing, interpolation, outcome-based selection,
statistical uncertainty bands, or trajectory exclusions. Export 183 x 120 mm.
"""
from __future__ import annotations
import argparse
import sys
from pathlib import Path
import numpy as np
import complete_task_evidence as E


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--qa-tools', type=Path, required=True)
    args = parser.parse_args()
    sys.path.insert(0, str(args.qa_tools))
    from audit_panel_alignment import require_matplotlib_panel_alignment
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    E.verify()
    plt.rcParams.update({'font.family': 'sans-serif',
                         'font.sans-serif': ['Microsoft YaHei', 'DejaVu Sans'],
                         'font.size': 8, 'axes.unicode_minus': False,
                         'axes.spines.top': False, 'axes.spines.right': False,
                         'axes.linewidth': .7, 'legend.frameon': False,
                         'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    methods = [('Best_fixed', '固定控制 B1', '#707070', '--'),
               ('EtR-A-Haninger', 'Haninger 规则适配', '#a88a53', ':'),
               ('Candidate1_v2', 'C1：模型优化', '#2878a5', '-'),
               ('Candidate2', 'C2：学习', '#ac557a', '-')]
    fig, axes = plt.subplots(2, 2, figsize=(7.2047, 4.7244))  # 183 x 120 mm
    fig.subplots_adjust(left=.105, right=.98, bottom=.12, top=.79,
                        wspace=.32, hspace=.55)
    ss = E.scenes('T2_motion')
    output = E.OUT / 'figures'; output.mkdir(exist_ok=True)
    for col, scene in enumerate((0, 2)):
        p, T = ss[scene]
        for name, label, color, style in methods:
            spec = next(s for s in E.specs('T2_motion') if s[0] == name and s[3] in (-1, 0))
            pol, update = E.build(spec); E.bind(pol, p)
            log, _, commands = E.episode(p, pol, T, E.SEED + scene, update)
            assert log is not None
            np.savez_compressed(output / f'push_scene{scene}_{name}.npz', log=log, commands=commands, p=p)
            axes[0, col].plot(log[:, 0], log[:, 3], label=label, color=color,
                              linestyle=style, linewidth=1.1)
            axes[1, col].plot(log[:, 0], (log[:, 2] - log[0, 2]) * 1000,
                              color=color, linestyle=style, linewidth=1.1)
        axes[0, col].axhline(4, color='#222222', linewidth=.6, linestyle='-.')
        material = '软' if col == 0 else '硬'
        axes[0, col].set_title(f'{chr(97 + col)}  均匀{material}组织：{p[0] / 1000:.0f} kPa',
                               loc='left', fontsize=8.5, pad=9)
        axes[1, col].set_title(f'{chr(99 + col)}  实际法向运动', loc='left', fontsize=8.5, pad=9)
        for row in range(2):
            ax = axes[row, col]
            ax.set_xlim(0, T)
            ax.set_xticks([0, 2, 4, 6, 8])
            ax.set_xlabel('时间（s）', labelpad=3)
            ax.tick_params(length=3, width=.6, labelsize=7)
    axes[0, 0].set_ylabel('接触力（N）')
    axes[1, 0].set_ylabel('相对初始位置的位移（mm）')
    axes[0, 1].set_ylabel('接触力（N）')
    axes[1, 1].set_ylabel('相对初始位置的位移（mm）')
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(.54, .98),
               ncol=2, fontsize=8, columnspacing=2., handlelength=2.4)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(fig, json_out=output / 'motion_response.alignment.json',
                                        tolerance_pt=1.5, gutter_tolerance_pt=1.5, strict=True)
    fig.savefig(output / 'motion_response.pdf')
    fig.savefig(output / 'motion_response.svg')
    fig.savefig(output / 'motion_response.png', dpi=300)
    fig.savefig(output / 'motion_response.tiff', dpi=600,
                pil_kwargs={'compression': 'tiff_lzw'})
    E.save('figure_metadata.json', dict(question='Material dependence of frozen-controller motion response',
             source_scenes=[0, 2], source_seed=E.SEED, neural_training_seed=0,
             excluded_methods='Four illustrative methods shown; all thirteen methods retained in primary tables',
             trajectory_exclusions=0, smoothing=False, axis_clipping=False,
             uncertainty='Individual deterministic scene/noise traces; no ensemble uncertainty implied',
             panel_roles=['soft force', 'stiff force', 'soft actual position', 'stiff actual position']))
    print(output / 'motion_response.png')


if __name__ == '__main__':
    main()
