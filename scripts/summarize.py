"""Run on the CPU allocation after Julia has written the source CSV tables."""
from pathlib import Path
import csv
import json
import platform
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter

root = Path(__file__).resolve().parents[1]
out = root / 'results'

def rows(name):
    with (out / name).open() as f:
        return list(csv.DictReader(f))

def slope(data, x, y):
    return float(np.polyfit(np.log([float(r[x]) for r in data]),
                            np.log([float(r[y]) for r in data]), 1)[0])

data = rows('convergence.csv')
summary = {'convergence': {}, 'kernel_coupling_orders': {},
           'python': platform.python_version(), 'numpy': np.__version__,
           'matplotlib': matplotlib.__version__}
styles = {'raw': ('#0072B2','o','-'), 'vertex': ('#D55E00','s','--'),
          'raw_richardson': ('#009E73','^','-.'),
          'vertex_richardson': ('#CC79A7','D',':')}
with plt.rc_context({'font.size': 10, 'axes.grid': True, 'grid.alpha': 0.2}):
    fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), layout='constrained')
    for ax, axis in zip(axes, ('real','imag','mixed')):
        summary['convergence'][axis] = {}
        for method, (color, marker, line) in styles.items():
            a = [r for r in data if r['axis']==axis and r['method']==method]
            if len(a)<2:
                continue
            a.sort(key=lambda r: float(r['step']))
            p=slope(a,'step','relative_error')
            tail=slope(a[:2],'step','relative_error')
            summary['convergence'][axis][method] = {
                'p_all':p,'p_finest_two':tail,
                'N':[int(r['N']) for r in a],
                'relative_errors':[float(r['relative_error']) for r in a]}
            ax.loglog([float(r['step']) for r in a],
                      [float(r['relative_error']) for r in a],
                      color=color,marker=marker,linestyle=line,
                      label=f'{method.replace("_", " ")} (p={p:.2f})')
        ax.set(title={'real':'Keldysh: factorized initial state',
                      'imag':'Imaginary time: equilibrium',
                      'mixed':'Mixed contour: correlated initial state'}[axis],
               xlabel='step (real or imaginary time)',ylabel='relative L2 Green function error')
        ticks=sorted({float(r['step']) for r in data if r['axis']==axis and r['method']=='raw'})
        ax.set_xticks(ticks,[f'{v:g}' for v in ticks])
        ax.xaxis.set_minor_formatter(NullFormatter())
        ax.legend(fontsize=7)
    for suffix in ('png','pdf'):
        fig.savefig(out / f'convergence.{suffix}',dpi=180)
    plt.close(fig)

kd=rows('kernel_coupling_scaling.csv')
for axis in ('real','imag'):
    for method in ('raw','vertex'):
        a=[r for r in kd if r['axis']==axis and r['method']==method]
        summary['kernel_coupling_orders'][f'{axis}_{method}']=slope(a,'V','absolute_kernel_error')

if (out/'link_if.csv').exists():
    summary['link_if']={}
    ld=rows('link_if.csv')
    for ed in sorted({r['ed'] for r in ld}):
        a=[r for r in ld if r['ed']==ed]
        if len(a)>=2:
            summary['link_if'][ed]={
                'order_in_N':-slope(a,'N','relative_error'),
                'max_split_mismatch':max(float(r['split_oracle_error']) for r in a),
                'max_trace_error':max(abs(complex(float(r['Z_re']),float(r['Z_im']))-1) for r in a)}
if (out/'link_if_strang.csv').exists():
    a=rows('link_if_strang.csv')
    summary['strang_link_if']={
        'order_in_step':-slope(a,'N','relative_error'),
        'max_split_mismatch':max(float(r['split_oracle_error']) for r in a),
        'max_trace_error':max(abs(complex(float(r['Z_re']),float(r['Z_im']))-1) for r in a),
        'relative_errors':[float(r['relative_error']) for r in a]}
if (out/'tti_kernel_audit.csv').exists():
    summary['tti_kernel_audit']=rows('tti_kernel_audit.csv')
if (out/'repaired_tti.csv').exists():
    summary['repaired_tti']=rows('repaired_tti.csv')
summary['figure_provenance']={
    'source':'convergence.csv',
    'error':'norm(g-g_reference)/norm(g_reference), all node samples included',
    'richardson':'2*g_fine_on_coarse_nodes - g_coarse, independently normalized',
    'uncertainty':'deterministic computation; no statistical error bars',
    'log_policy':'all plotted source errors are strictly positive',
    'alt_text':'Three panels compare raw, shifted, and Richardson-extrapolated Green functions. The mixed shifted scheme approaches a nonzero error; other trends are in summary.json.',
    'scope':'Single spinless level, one bath mode, beta=1, impurity energy=0.7, bath energy=1, squared coupling=0.5.'}
(out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(summary,ensure_ascii=False,indent=2))
