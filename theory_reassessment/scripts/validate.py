"""Bounded result assertions and reproducible summaries; execute on CPU server."""
from foundation import *
from contour_audit import wrong_continuum, source_if


def read(name):
    return list(csv.DictReader(open(ROOT/'results'/name)))


def complexcol(row, prefix):
    return float(row[prefix+'_re'])+1j*float(row[prefix+'_im'])


def main():
    import scipy
    out={'environment':dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__)}
    cases={c['id']:c for c in read('cases.csv')}
    blocks=read('source_blocks.csv')
    out['source_blocks_max_absolute_error']=max(float(r['error']) for r in blocks)
    assert out['source_blocks_max_absolute_error']<1e-12
    out['source_blocks_checked']=len(blocks)
    algebra=read('local_algebra.csv')
    out['square_identity_max_error']=max(float(r['error']) for r in algebra if r['test']=='square')
    out['local_orders']=[]
    for r in ('1','2'):
        for shared in ('False','True'):
            for method in ('q','normal','link','weyl'):
                a=sorted((x for x in algebra if x['r']==r and x['shared']==shared and x['test']==method),key=lambda x:float(x['h']))
                if len(a)<2:continue
                p=np.log(float(a[1]['error'])/float(a[0]['error']))/np.log(2)
                out['local_orders'].append(dict(r=r,shared=shared,method=method,p=float(p)))
                assert p>(2.9 if method=='weyl' else 1.9)
    conv=read('operator_convergence.csv');out['global_orders']=[]
    for r in ('1','2'):
        for axis in ('real','imag','mixed'):
            for method in ('q_noextra','link','weyl'):
                a=sorted((x for x in conv if x['r']==r and x['axis']==axis and x['method']==method),key=lambda x:int(x['N']))
                p=np.log(float(a[-2]['error'])/float(a[-1]['error']))/np.log(2)
                out['global_orders'].append(dict(r=r,axis=axis,method=method,p=float(p),error128=float(a[-1]['error'])))
                if method=='weyl':assert 1.85<p<2.15
    out['occupation_leading']=[r for r in read('leading_error.csv') if r['N']=='512']
    links=read('link_and_sign.csv');gauss=read('gaussian.csv')
    out['main_contour_errors']=[r for r in links if r['id'].startswith('main') and int(r['N']) in (16,32,64,128)]
    out['main_native_errors']=[r for r in gauss if r['id'].startswith('main') and int(r['N']) in (16,32,64,128)]
    case=dict(cases['main_mixed_128_128']);case['N']='4096'
    wrong,_=wrong_continuum(case);ref=physical_reference(case)
    out['mixed_wrong_continuum_relative_error']=float(norm(wrong-ref)/norm(ref))
    # Verify original greater Green functions with full 4x4 finite operators.
    leading=read('green_leading.csv');err=0.
    for N in (16,32,64,128):
        case=cases[f'main_real_{N}_{N}'];nodes,ds=load_case(case)
        a=[r for r in leading if r['model']=='main' and int(r['N'])==N]
        for mode,prefix in (('raw','raw'),('vertex','link')):
            g,_,_=gaussian(case,nodes,ds[mode]);gg=np.array([complexcol(r,prefix) for r in a])
            err=max(err,float(norm(g-gg)))
    out['green_operator_vs_gaussian_max_l2']=err
    assert err<1e-10
    a=[r for r in leading if r['model']=='main' and r['N']=='128']
    ref=np.array([complexcol(r,'reference') for r in a]);den=norm(ref)
    out['green_leading_norms']={key:float(norm([complexcol(r,p) for r in a])/den) for key,p in (('Q_bulk','qbulk'),('endpoint','endpoint'),('link','linkcoef'))}
    out['green_leading_norms']['raw_total']=float(norm([complexcol(r,'qbulk')+complexcol(r,'endpoint') for r in a])/den)
    out['green_leading_residual128']={p:max(float(r[p+'_coefficient_residual']) for r in a) for p in ('raw','link')}
    # Native package contraction against independently exponentiated matrices.
    native=read('native.csv');eg=ez=0.
    for axis in ('real','imag','mixed'):
        for N in (2,4,8,16):
            case=cases[f'main_{axis}_{N}_{N}'];nodes,ds=load_case(case)
            for mode in ('raw','vertex'):
                g,Z,_=gaussian(case,nodes,ds[mode])
                rr=[r for r in native if r['axis']==axis and int(r['N'])==N and r['method']==mode]
                assert len(rr)==len(g),(axis,N,mode,len(rr),len(g))
                eg=max(eg,max(abs(complexcol(r,'g')-g[int(r['k'])-1]) for r in rr))
                ez=max(ez,max(abs(complexcol(r,'Z')-Z) for r in rr))
    out['native_max_g_error']=float(eg);out['native_max_Z_error']=float(ez)
    assert max(eg,ez)<1e-9
    pkg=read('package.csv');assert len(pkg)==34
    out['package_cases']=len(pkg)
    out['package_same_target_max_error']=max(float(r['same_target_error']) for r in pkg)
    out['package_K_coefficient_max_error']=max(float(r['K_coefficient_error']) for r in pkg if r['K_coefficient_error']!='NaN')
    assert out['package_same_target_max_error']<1e-9
    out['package_errors']=pkg
    out['package_orders']=[]
    for axis in ('real','imag','mixed'):
        for method in ('q','weyl'):
            a=sorted((r for r in pkg if r['r']=='1' and r['green']=='false' and r['axis']==axis and r['method']==method),key=lambda r:int(r['N']))
            p=float(np.log(float(a[-2]['physical_error'])/float(a[-1]['physical_error']))/np.log(2))
            out['package_orders'].append(dict(axis=axis,method=method,p=p,error16=float(a[-1]['physical_error'])))
            if method=='weyl':assert 1.95<p<2.05
    alt=read('alternatives.csv');assert len(alt)==12
    out['alternative_native_cases']=len(alt)
    out['alternative_native_max_error']=max(float(r[k]) for r in alt for k in ('g_error','Z_error'))
    assert out['alternative_native_max_error']<1e-9
    out['warm_timings']=read('warm_timings.csv')
    half=[r for r in out['occupation_leading'] if r['method']=='half_endpoints']
    assert len(half)==3
    out['half_occupation_leading_residual512']=max(abs(float(r['error_over_h'])-float(r['predicted'])) for r in half)
    assert out['half_occupation_leading_residual512']<1e-4
    with open(ROOT/'results/validation.json','w') as f:json.dump(out,f,indent=2,ensure_ascii=False)
    print(json.dumps({k:v for k,v in out.items() if k not in ('main_contour_errors','main_native_errors','package_errors','occupation_leading')},indent=2))


if __name__=='__main__':main()
