"""picks (T ms, Vavg) -> Vint (formula de velocidade media em tempo) -> modelo v(z) em ONL 25 m"""
import sys, json, numpy as np
from mkmod import wonl
import pick
def vint_from_picks(P, twb=1333.33, vw=1500.):
    P=[(twb,vw)]+[(t,v) for t,v in P if t>twb+50]
    rows=[]
    for (t1,v1),(t2,v2) in zip(P[:-1],P[1:]):
        T1,T2=t1/1000,t2/1000; vi=(v2*T2-v1*T1)/(T2-T1)
        rows.append(dict(T1=t1,T2=t2,V1=v1,V2=v2,z1=v1*T1/2,z2=v2*T2/2,Vint=vi))
    return rows
def model_from_rows(rows, nx=512, nz=200, dz=25., vw=1500.):
    z=np.arange(nz)*dz; v=np.full(nz,vw)
    for r in rows: v[z>=r['z1']]=r['Vint']
    return np.tile(v,(nx,1))
if __name__=='__main__':
    vfix, npt, out = sys.argv[1], int(sys.argv[2]), sys.argv[3]
    ev=[(t,v) for t,v,n in pick.events(pick.picks(vfix,npt,tmin=1600.)) if t<3900]
    rows=vint_from_picks(ev)
    wonl(out, model_from_rows(rows), 25.)
    json.dump(dict(picks=ev,rows=rows), open(out+'.json','w'), indent=1)
    for r in rows: print('T %5.0f-%5.0f ms  Vavg %5.0f->%5.0f  z %5.0f-%5.0f m  Vint %5.0f'%(r['T1'],r['T2'],r['V1'],r['V2'],r['z1'],r['z2'],r['Vint']))
