"""Explicit native-lattice targets for consistent links and half endpoints."""
from foundation import *
from contour_audit import source_if
from package_inputs import jl,jmat


def run():
    out=[]
    for c in csv.DictReader(open(ROOT/'results/cases.csv')):
        if not c['id'].startswith('main') or c['axis']=='imag' or int(c['N']) not in (4,8,16):continue
        nodes,ds=load_case(c)
        for method in ('derived_link','half_endpoints'):
            D,_=source_if(c,nodes,links=method=='derived_link',half=method=='half_endpoints')
            g,Z,_=gaussian(c,nodes,D)
            ns='['+','.join(f'(:{b},{j})' for b,j in nodes)+']'
            out.append('(;axis=:'+c['axis']+',N='+c['N']+',method=:'+method+',nodes='+ns+',D='+jmat(D)+',g=['+','.join(jl(z) for z in g)+'],Z='+jl(Z)+')')
    (ROOT/'results/alternative_cases.jl').write_text('const CASES = [\n'+',\n'.join(out)+'\n]\n')
    print('generated',len(out),'native alternative targets')


if __name__=='__main__':run()
